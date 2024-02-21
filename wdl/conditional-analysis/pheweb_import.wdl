version 1.0

workflow conditional_pheweb {

  input {
    String docker
    String release
  }
  call split_inputs {input:docker=docker}

  scatter (i in range(length(split_inputs.chain_chunks))) {
    call pheweb_import_munge {
      input :
      docker = docker,
      release = release,
      chains = split_inputs.chain_chunks[i],
      outputs = split_inputs.output_chunks[i]
    }
  }

  call merge_pheweb {
    input :
    docker = docker,
    release = release,
    sql_chunks = pheweb_import_munge.csv_sql,
    output_chunks = flatten(pheweb_import_munge.munged_regenie),
  }
  
}


task merge_pheweb {
  input {
    Array[File] sql_chunks
    Array[String] output_chunks
    String docker
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
    cpu: "4"
    docker: "${docker}"
    memory: "4 GB"
    disks: "local-disk ~{disk_size} HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: "1"
    
  }
}

task pheweb_import_munge{
  input {
    String release
    File regions
    File chains
    File outputs
    String docker
    
  }
  Array[File] cond_locus_hits = read_lines(chains)
  Array[File] regenie_outputs = read_lines(outputs)

  Int disk_size = ceil(size(cond_locus_hits[0],"MB")*length(cond_locus_hits)*2/1000 + size(regenie_outputs[0],"MB")*length(regenie_outputs)*2/1000 +10)
  String prefix = "finngen_R" + release
  String out_file = prefix + "_sql.merged.txt"

  command <<<
    python3 <<CODE
    import os,sys,math
    
    release = '~{release}'
    prefix = '~{prefix}'
    out_dir = "."
    region_file = '~{regions}'
    hits = '~{write_lines(cond_locus_hits)}'
  
    out_file = os.path.join(out_dir,f"{prefix}_sql.merged.txt")
    #reads in all paths
    with open(hits) as f:hits_paths = [elem.strip() for elem in f.readlines()]
    #loop over region data, find matching file(s), merge info and build sql table
    with open(region_file) as f,open(out_file,'wt') as o:
        count = 0
        for line in f:
            pheno,chrom,region,locus,*_ = line.strip().split()
            id_string = f"{prefix}_{pheno}_{locus}"
            matches =[path for path in hits_paths if id_string in path]
            if matches:
                assert len(matches)==1
                count +=1
                sys.stdout.write("\r%s" % f"{pheno} {count}/{len(hits_paths)}                                                  ")
                sys.stdout.flush()
                hit_file = matches[0]
                file_root = os.path.basename(hit_file).split("chr")[0]
                start,end=region.split(':')[1].split("-")
                with open(hit_file) as i: variants = [elem.strip().split()[0] for elem in i][1:]
                out_data = [release,"conditional",pheno,chrom,start,end,len(variants),"",','.join(variants),file_root]
                o.write(','.join([f'"{elem}"' for elem in out_data])+'\n')


    #fixing columns of regenie outputs
    inputs = '~{write_lines(regenie_outputs)}'
    separator = " "
    input_columns = ['CHROM', 'GENPOS', 'ID', 'ALLELE0', 'ALLELE1', 'A1FREQ', 'INFO', 'N', 'TEST', 'BETA', 'SE', 'CHISQ', 'LOG10P', 'EXTRA']
    out_columns = ['SNPID','CHR','rsid','POS','Allele1','Allele2','AF_Allele2','p.value_cond','BETA_cond','SE_cond']
    map_columns = {"SNPID":"ID","CHR":"CHROM","rsid":"ID","POS":"GENPOS","Allele1":"ALLELE0","Allele2":"ALLELE1","AF_Allele2":"A1FREQ","p.value_cond":"LOG10P","BETA_cond":'BETA','SE_cond':'SE'}
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
    ls && echo ~{disk_size}
  >>>
  output {
    File csv_sql = out_file
    Array[File] munged_regenie = glob("./finngen*conditional")    
  }
  
  runtime {
    cpu: "4"
    docker: "${docker}"
    memory: "4 GB"
    disks: "local-disk ~{disk_size} HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: "1"
      
  }
  

}

task split_inputs {
  input {
    String docker
    File chains
    File outputs
    Int chunks
  }

  command <<<
  split ~{chains} -dn l/~{chunks} chains
  split ~{outputs} -dn l/~{chunks} outputs
  >>>

  output {
    Array[String] chain_chunks = glob("./chains*")
    Array[String] output_chunks = glob("./outputs*")
  }
  
  runtime {
     cpu: "4"
     docker: "${docker}"
     memory: "4 GB"
     disks: "local-disk 2 HDD"
     zones: "europe-west1-b europe-west1-c europe-west1-d"
     preemptible: "1"
      
  }
}