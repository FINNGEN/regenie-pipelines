version 1.0

workflow conditional_pheweb {

  input {
    String release
    File chains
    File outputs
    File regions
    Int chunks
  }

  call split_inputs {input: chains=chains, outputs=outputs, chunks=chunks}

  # chains is chunked (not regions) because each chunk's chain paths get localized as files --
  # chunking regions instead and passing the full chains list to every shard would localize every
  # chain file once per shard. regions stays whole; build_sql_chunk subsets it down to just the
  # loci present in its chain chunk.
  scatter (chain_chunk in split_inputs.chain_chunks) {
    call build_sql_chunk {
      input :
      release = release,
      regions = regions,
      chains = chain_chunk,
      all_outputs = outputs
    }
  }

  scatter (output_chunk in split_inputs.output_chunks) {
    call munge_outputs_chunk {input: outputs=output_chunk}
  }

  call merge_pheweb {
    input :
    release = release,
    sql_chunks = build_sql_chunk.csv_sql,
    output_chunks = flatten(munge_outputs_chunk.munged_regenie),
  }

}


task split_inputs {
  input {
    File chains
    File outputs
    Int chunks
  }

  command <<<
    set -euo pipefail
    split ~{chains} -den l/~{chunks} chains_
    split ~{outputs} -den l/~{chunks} outputs_
  >>>

  output {
    Array[File] chain_chunks = glob("./chains_*")
    Array[File] output_chunks = glob("./outputs_*")
  }

  runtime {
    disks: "local-disk 2 HDD"
  }
}

task merge_pheweb {
  input {
    Array[File] sql_chunks
    Array[String] output_chunks
    String release
  }

  String out_file = "R" + release +"_sql.merged.txt"
  String outputs_list = "R" + release +"_sql.outputs.txt"
  Int disk_size = ceil(size(sql_chunks[0],"MB")*length(sql_chunks)*2/1000 + size(output_chunks[0],"MB")*length(output_chunks)*2/1000 +10)
  command <<<
  echo ~{disk_size}
  while read f;
  do cat $f >> ~{out_file}
  done < ~{write_lines(sql_chunks)}

  cat ~{write_lines(output_chunks)} > ~{outputs_list}
  >>>

  output {
    File sql = out_file
    File outputs = outputs_list
  }

  runtime {
    disks: "local-disk ~{disk_size} HDD"
  }
}

task build_sql_chunk {
  input {
    String release
    File regions
    File chains
    File all_outputs
  }
  Array[File] cond_locus_hits = read_lines(chains)

  Int disk_size = ceil(size(cond_locus_hits[0],"MB")*length(cond_locus_hits)*2/1000 +10)
  String prefix = "finngen_R" + release
  String out_file = prefix + "_sql.merged.txt"

  command <<<
    python3 <<CODE
    import os,sys

    release = '~{release}'
    prefix = '~{prefix}'
    out_dir = "."
    region_file = '~{regions}'
    hits = '~{write_lines(cond_locus_hits)}'
    all_outputs_file = '~{all_outputs}'

    out_file = os.path.join(out_dir,f"{prefix}_sql.merged.txt")
    #reads in all paths
    with open(hits) as f:hits_paths = [elem.strip() for elem in f.readlines()]
    with open(all_outputs_file) as f:conditional_paths = [elem.strip() for elem in f.readlines()]

    # this chunk's chain files were already localized (chains is chunked, not regions) -- filter
    # the full region_file down to just this chunk's loci in one linear pass, instead of scanning
    # every region against hits_paths (quadratic overall). Match on the "{pheno}_{locus}" suffix
    # of the chain filename rather than reconstructing the full "{prefix}_{pheno}_{locus}" id --
    # the release-derived prefix isn't guaranteed to match how upstream actually named the file
    # (e.g. special releases), but the pheno/locus suffix always does.
    chunk_ids = set()
    for path in hits_paths:
        id_string = os.path.basename(path)
        if id_string.endswith(".independent.snps"):
            id_string = id_string[:-len(".independent.snps")]
        chunk_ids.add(id_string)

    relevant_regions = []
    matched_ids = set()
    with open(region_file) as f:
        for line in f:
            pheno,chrom,region,locus,*_ = line.strip().split()
            suffix = f"_{pheno}_{locus}"
            matches = [id_string for id_string in chunk_ids if id_string.endswith(suffix)]
            if matches:
                assert len(matches)==1, f"ambiguous chain match for {suffix}: {matches}"
                id_string = matches[0]
                relevant_regions.append((pheno,chrom,region,locus,id_string))
                matched_ids.add(id_string)

    # chains/regions are always produced together upstream -- a chain file in this chunk with no
    # matching region entry is a real pipeline inconsistency, fail loudly rather than skip it.
    missing = chunk_ids - matched_ids
    assert not missing, f"chain file(s) with no matching region entry: {missing}"

    with open(out_file,'wt') as o:
        count = 0
        for pheno,chrom,region,locus,id_string in relevant_regions:
            matches = [path for path in hits_paths if id_string in path]
            assert len(matches)==1
            hit_file = matches[0]
            count +=1
            sys.stdout.write("\r%s" % f"{pheno} {count}/{len(hits_paths)}                                                  ")
            sys.stdout.flush()
            file_root = os.path.basename(hit_file).split("chr")[0]
            start,end=region.split(':')[1].split("-")
            with open(hit_file) as i: variants = [elem.strip().split()[0] for elem in i][1:]

            # a naturally-converged chain has one more step file than hit rows (the final,
            # non-significant confirming step); a chain truncated by MAX_STEPS doesn't -- drop
            # the last, unconfirmed hit in that case rather than claim it as terminal.
            n_conditional_files = len([p for p in conditional_paths if id_string in p])
            if n_conditional_files < len(variants):
                variants = variants[:-1]

            out_data = [release,"conditional",pheno,chrom,start,end,len(variants),"",','.join(variants),file_root]
            o.write(','.join([f'"{elem}"' for elem in out_data])+'\n')
    CODE
  >>>

  output {
    File csv_sql = out_file
  }

  runtime {
    disks: "local-disk ~{disk_size} HDD"
  }
}

task munge_outputs_chunk {
  input {
    File outputs
  }
  Array[File] regenie_outputs = read_lines(outputs)

  Int disk_size = ceil(size(regenie_outputs[0],"MB")*length(regenie_outputs)*2/1000 +10)

  command <<<
    python3 <<CODE
    import os,sys,math

    out_dir = "."
    #fixing columns of regenie outputs
    separator = " "
    input_columns = ['CHROM', 'GENPOS', 'ID', 'ALLELE0', 'ALLELE1', 'A1FREQ', 'INFO', 'N', 'TEST', 'BETA', 'SE', 'CHISQ', 'LOG10P', 'EXTRA']
    out_columns = ['SNPID','CHR','rsid','POS','Allele1','Allele2','AF_Allele2','p.value_cond','BETA_cond','SE_cond']
    map_columns = {"SNPID":"ID","CHR":"CHROM","rsid":"ID","POS":"GENPOS","Allele1":"ALLELE0","Allele2":"ALLELE1","AF_Allele2":"A1FREQ","p.value_cond":"LOG10P","BETA_cond":'BETA','SE_cond':'SE'}
    inputs = '~{write_lines(regenie_outputs)}'
    with open(inputs) as f:paths = [elem.strip() for elem in f.readlines()]
    n_paths = len(paths)
    for i,path in enumerate(paths):
        sys.stdout.write("\r%s" % f"{path} {i+1}/{n_paths}                                                  ")
        sys.stdout.flush()
        out_file = os.path.join(out_dir,os.path.basename(path))
        with open(path) as i,open(out_file,'wt') as o:
            header_index = {h:i for i,h in enumerate(next(i).strip().split(separator))}
            o.write(separator.join(out_columns) +'\n')
            for line in i:
                line = line.strip().split(separator)
                out_line = []
                for key in out_columns:
                    #i take the column mapping and then i get the data from the input header mapping in return
                    data_index = header_index[map_columns[key]]
                    value = str(line[data_index])                    #get data from input
                    if key =="p.value_cond":value =  math.pow(10,-float(value)) #fix pval
                    out_line.append(value)

                o.write(separator.join(map(str,out_line)) +'\n')
    CODE
  >>>
  output {
    Array[File] munged_regenie = glob("./finngen*conditional")
  }

  runtime {
    disks: "local-disk ~{disk_size} HDD"
  }
}
