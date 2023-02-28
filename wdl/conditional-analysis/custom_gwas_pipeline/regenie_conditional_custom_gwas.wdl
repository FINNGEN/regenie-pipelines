version 1.0

workflow conditional_analysis {

    input {
        String docker
        String phenotype
        String release
        Boolean is_binary

        String mlogp_col
        String pval_col
        String chr_col
        String pos_col
        String ref_col
        String alt_col
        Float conditioning_mlogp_threshold
        Float locus_mlogp_threshold  
        Array[String] chroms
        String covariates
        File cov_file
        File pheno_file
        File sumstats
        File null
        Int threshold_cov_count

    }
    String proper_covars = sub(covariates,"PC\\{1:10\\}","PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10")
    String prefix = "CUSTOM_GWAS"

    
    call join_pheno_cov{
        input:pheno=phenotype,
        pheno_file=pheno_file,
        cov_file=cov_file,
        docker=docker
    }
    call region_selection{
        input: pheno=phenotype, 
        sumstats=sumstats, 
        chromosome_col=chr_col,
        position_col=pos_col,
        allele1_col=ref_col,
        allele2_col=alt_col,
        p_col=pval_col,
        freq_col="af_alt",
        se_col="sebeta",
        beta_col="beta"

    }
    if(region_selection.had_results){
        call extract_cond_regions {
        input: pheno=phenotype,
        region = region_selection.bed,
        mlogp_threshold = locus_mlogp_threshold,
        docker=docker,
        mlogp_col = mlogp_col,
        chr_col=chr_col,
        pos_col = pos_col,
        ref_col=ref_col,
        alt_col=alt_col,
        chroms=chroms,
        sumstats=sumstats
    }



    Array[Array[String]] all_regions = read_tsv(extract_cond_regions.gw_sig_res)

    # loop over all regions running one variant per shard
    scatter (region in all_regions) {
        String npheno = region[0]
        String chrom = region[1]
        String locus_region = region[2] + " " + region[3]
    call regenie_conditional {
        input: docker = docker,
        prefix=prefix,
        locus_region=locus_region,
        pheno=npheno,
        chrom=chrom,
        covariates = proper_covars,
        mlogp_col = mlogp_col,
        chr_col=chr_col,
        pos_col = pos_col,
        ref_col=ref_col,
        alt_col=alt_col,
        pval_threshold=conditioning_mlogp_threshold,
        sumstats=sumstats,
        pheno_file=join_pheno_cov.joined_cov_pheno,
        null=null
    }    
    }

    Array[File] results = flatten(regenie_conditional.conditional_chains)
    #There is no reason for a complicated task, since what it does is join the pheno-specific regenion files together.
    # I am joining all files together, so it's the same 
    call cg_merge_results{input:docker=docker,files = results,phenotype=phenotype}

    #merges all regions into single file
    Array[File] regenie_outputs = flatten(regenie_conditional.regenie_output)
    call pheweb_import_munge{
        input:docker=docker,
        prefix=prefix,
        release=release,
        pheno=phenotype,
        cond_locus_hits=results,
        regions=extract_cond_regions.gw_sig_res,
        regenie_outputs=regenie_outputs
    }
    }
    

    output{
        String sql_import_file = select_first([pheweb_import_munge.csv_sql,"None"])
        String independent_snps = select_first([cg_merge_results.pheno_independent_snps,"None"])
        Array[String] conditional_data = select_first([regenie_outputs,[]])
        Boolean had_results = region_selection.had_results
    }
}

task join_pheno_cov{
    input{
        String pheno
        File pheno_file
        File cov_file
        String docker
    }
    command <<<
    set -eux
    #for phenotype file, only take FID and phenotype column
    pheno_col=$(zcat -f ~{pheno_file}|head -n1|tr "\t" "\n"|sed -n "/^~{pheno}$/ =")
    join --header -1 1 -2 1 -t $'\t' \
        <(cat <(zcat -f ~{cov_file}|head -n1) <(zcat -f ~{cov_file}|tail -n+2|sort -k1)) \
        <(cat <(zcat -f ~{pheno_file}|head -n1) <(zcat -f ~{pheno_file}|tail -n+2|sort -k1)|cut -f 1,$pheno_col) \
        |gzip > joined_cov_pheno.gz
    >>>

    output{
        File joined_cov_pheno = "joined_cov_pheno.gz"

    }

    runtime{
        docker: "~{docker}"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 10 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}

task region_selection {
    #NOTE: regular release region selection
    input {
        String pheno
        File sumstats
        String docker
        String chromosome_col
        String position_col
        String allele1_col
        String allele2_col
        String beta_col
        String se_col
        String p_col
        String freq_col
        Int window
        Int max_region_width
        Float window_shrink_ratio
        Float p_threshold
        Float? minimum_pval
        String? set_variant_id_map_chr
    }
    Boolean scale_se_by_pval = false
    Boolean x_chromosome = true
    Boolean set_variant_id = true
    String rsid_col = ""
    
    String delimiter = "TAB"
    
    # can be helpful if adding finemapping with relaxed threshold after more stringent has already ben run.
    # does not include regions with lead snp < this
    

    command <<<

        make_finemap_inputs.py \
            --sumstats ~{sumstats} \
            --rsid-col "~{rsid_col}" \
            --chromosome-col ~{chromosome_col} \
            --position-col "~{position_col}" \
            --allele1-col "~{allele1_col}" \
            --allele2-col "~{allele2_col}" \
            --freq-col "~{freq_col}" \
            --beta-col "~{beta_col}" \
            --se-col "~{se_col}" \
            --p-col "~{p_col}" \
            --delimiter "~{delimiter}" \
            --grch38 \
            --exclude-MHC \
            --no-upload \
            --prefix ~{pheno} \
            --out ~{pheno} \
            --window ~{window} \
            --max-region-width ~{max_region_width} \
            --window-shrink-ratio ~{window_shrink_ratio} \
            ~{true='--scale-se-by-pval ' false=' ' scale_se_by_pval} \
            ~{true='--x-chromosome' false=' ' x_chromosome} \
            ~{true='--set-variant-id ' false=' ' set_variant_id} \
            ~{true='--set-variant-id-map-chr ' false=' ' defined(set_variant_id_map_chr)}~{set_variant_id_map_chr} \
            --p-threshold ~{p_threshold} \
            ~{true='--min-p-threshold ' false='' defined(minimum_pval)}~{minimum_pval} \
            --wdl

            res=`cat ~{pheno}_had_results`

            if [ "$res" == "False" ]; then
                touch ~{pheno}".z"
                touch ~{pheno}".lead_snps.txt"
                touch ~{pheno}".bed"
            fi
    >>>

    output {

        Array[File] zfiles = glob("*.z")
        File bed = pheno + ".bed"
        File log = pheno + ".log"
        Boolean had_results = read_boolean("~{pheno}_had_results")
    }

    runtime {
        docker: "~{docker}"
        cpu: "1"
        memory: "60 GB"
        disks: "local-disk 20 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}



task pheweb_import_munge{
  input {
    String release
    String prefix
    String pheno
    File regions
    Array[File] cond_locus_hits
    String docker
    Array[File] regenie_outputs
  }
  Int disk_size = ceil(size(cond_locus_hits,'GB')) + ceil(size(regenie_outputs,'GB')) + 10
  
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
    ls
  >>>
  output {
    File csv_sql = out_file
    Array[File] munged_regenie = glob("./finngen*conditional")
  }
  
  runtime {
    cpu: "4"
    docker: "${docker}"
    memory: "4 GB"
    disks: "local-disk ${disk_size} HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: "1"
      
  }
  

}

task cg_merge_results {
    input{
        Array[File] files
        String docker
        String phenotype
    }

    command <<<
        #write files 
        head -n1 ~{files[0]} >> ~{phenotype}.independent_snps
        while read f; do tail -n+2 $f >> ~{phenotype}.independent_snps.txt; done <~{write_lines(files)}
    >>>

    output{
        File pheno_independent_snps = phenotype+".independent_snps.txt"
    }

    runtime {
    cpu: "1"
    docker: "${docker}"
    memory: "4 GB"
    disks: "local-disk 10 HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: "1"
  
  }
}

task merge_results {

  input {
    Array[File] result_list
    File phenos_list
    String docker
    String prefix
  }

  command <<<
    cat ~{write_lines(result_list)} > results.txt 
    # add prefix and suffix to make sure matching is not based on partially similar pheno names (ICD roots)
    cat ~{phenos_list} |awk  '{print "'~{prefix}'_"$0"_chr"}' > phenos.txt
    
    # match pheno names and files and extract phenos back by removing pre/suff ix
    grep -of phenos.txt results.txt  | sed 's/_chr//g' > match_phenos.txt 
    # same grepping but only keeping file path
    grep -f phenos.txt results.txt > match_paths.txt 

    #now paste together the two files to have a tab separated pheno and filepath
    paste match_phenos.txt match_paths.txt > tmp.txt && head tmp.txt
    
    # write header in each pheno file
    head -n1 results.txt | xargs head -n1 > header.txt
    # loop through phenos and write header line to each pheno output
    while read f; do cp header.txt $f.independent_snps.txt; done < <(cut -f1 tmp.txt)
    # append to pheno output  the content of each matching regenie output
    while read line; do arr=($line) && echo ${arr[0]} && cat ${arr[1]} |sed -e 1d >> ${arr[0]}.independent_snps.txt; done < tmp.txt
  >>>
  
  output {
    Array[File] pheno_independent_snps = glob("${prefix}*.txt")    
  }
  
  runtime {
    cpu: "1"
    docker: "${docker}"
    memory: "4 GB"
    disks: "local-disk 10 HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: "1"
  
  }
   
}


task regenie_conditional {
  
  input {
    File sumstats
    # GENERAL PARAMS
    String docker
    String prefix
    # hit info
    String locus_region
    String pheno
    String chrom
    # files to localize 
    File pheno_file
    String bgen_root
    File null
    File sumstats
    # column names and stuff
    String chr_col
    String pos_col
    String ref_col
    String alt_col
    String mlogp_col
    String beta
    String sebeta
    # Script parameters/options
    Float pval_threshold
    Int max_steps
    String covariates
    String? regenie_params
    Int cpus
  }

  # localize all files based on roots and pheno/chrom info
  
  File sum_tabix = sumstats + ".tbi"
  
  File bgen = sub(bgen_root,'CHROM',chrom)
  File bgen_sample = bgen + ".sample"
  File bgen_index = bgen + ".bgi"


  # runtime params based on file sizes
  Int disk_size = ceil(size(bgen,'GB')) + ceil(size(sumstats,'GB')) + ceil(size(null,'GB')) + ceil(size(pheno_file,'GB')) + 1

  command <<<
    
    echo ~{pheno} ~{chrom} ~{cpus} 
    tabix -h ~{sumstats}  ~{chrom} > region_sumstats.txt

    python3 /scripts/regenie_conditional.py \
    --out ./~{prefix}  \
    --bgen ~{bgen}  \
    --null-file ~{null}  \
    --sumstats region_sumstats.txt \
    --pheno-file ~{pheno_file} \
    --pheno ~{pheno} \
    --locus-region ~{locus_region}  \
    --pval-threshold ~{pval_threshold} \
    --max-steps ~{max_steps} \
    --chr-col ~{chr_col} \
    --pos-col ~{pos_col} \
    --ref-col ~{ref_col} \
    --alt-col ~{alt_col} \
    --mlogp-col ~{mlogp_col} \
    --beta-col ~{beta} \
    --sebeta-col ~{sebeta} \
    --covariates ~{covariates} \
    ~{if defined(regenie_params) then " --regenie-params " + regenie_params else ""} \
    --log info

  >>>
  output {
    Array[File] conditional_chains = glob("./${prefix}*.snps")
    Array[File] logs = glob("./${prefix}*.log")
    Array[File] regenie_output = glob("./${prefix}*.conditional")    
  }
  
  runtime {
    cpu: "~{cpus}"
    docker: "${docker}"
    memory: "4 GB"
    disks: "local-disk ${disk_size} HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: "1"
  
  }
  
}


task filter_covariates {

  input {
    File pheno_file
    Array[String] covariates
    File pheno_list
    String docker
    Int threshold_cov_count
    }

    String outfile = "./pheno_cov_map_" + threshold_cov_count + ".txt"
    Int disk_size = ceil(size(pheno_file,'GB')) + 2 * 2
    
    command <<<

      set -euxo pipefail
      
      python3 <<CODE
      
      import pandas as pd
      import numpy as np
      
      #read in phenos as list of phenos regardless
      tot_phenos = []
      phenos_groups = []
      with open('~{pheno_list}') as i:
          for line in i:
              phenos = line.strip().split()
              phenos_groups.append(phenos)
              tot_phenos += phenos    

      #read in phenos mapping all valid entries to 1 and NAs to 0
      pheno_df= pd.read_csv('~{pheno_file}',sep='\t',usecols=tot_phenos).notna().astype(int)
      print(pheno_df)
      # read in covariates getting absolute values (handles PCs)
      covariates= '~{sep="," covariates}'.split(',')
      cov_df= pd.read_csv('~{pheno_file}',sep='\t',usecols=covariates).abs()
      print(cov_df)

      # now for each pheno calculate product of each covariate with itself
      
      with open('~{outfile}','wt') as o,open('~{outfile}'.replace('.txt','.err.txt'),'wt') as tmp_err:
          for i,pheno_list in enumerate(phenos_groups):
              pheno_name = ','.join(pheno_list)
              # for each group of phenos (possibly a single one) multiply all covs and pheno column and count how many non 0 entries are there: this means that the entry has a valid pheno and a non null covariates.
              df = pd.DataFrame()
              for pheno in pheno_list:
                  m = (cov_df.mul(pheno_df[pheno],0)>0).sum().to_frame(pheno)
                  df = pd.concat([df,m],axis =1)

              print(f"{i+1}/{len(phenos_groups)} {pheno_name}")
              #If it's a group of phenos the min function will return the lowest count across all phenos
              tmp_df = df[pheno_list].min(axis =1)
              covs = tmp_df.index[tmp_df >= ~{threshold_cov_count}].tolist()
              missing_covs = [elem for elem in covariates if elem not in covs]
              if missing_covs:tmp_err.write(f"{pheno_name}\t{','.join(missing_covs)}\n")
              o.write(f"{pheno_name}\t{','.join(covs)}\n")
      
      CODE

    >>>
      output {
	File cov_pheno_map = outfile
	File cov_pheno_warning = sub(outfile,'.txt','.err.txt')
      }
  
  runtime {
    cpu: "1"
    docker: "${docker}"
    memory: "64 GB"
    disks: "local-disk ${disk_size} HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: "1"
  }

}

  
task extract_cond_regions {
  
  input {
    
    String pheno
    File sumstats
    File region
    String chr_col
    String pos_col
    String ref_col
    String alt_col
    String mlogp_col
    Array[String] chroms 

    String docker
    Float mlogp_threshold

  }
  
  String outfile= pheno + "_sig_hits.txt" 
  
  command <<<
    #get column numbers
    chrnum=$(zcat ~{sumstats} |head -n1|tr "\t" "\n"|cat -n|grep ~{chr_col}|cut -f1|tr -d " ")
    posnum=$(zcat ~{sumstats} |head -n1|tr "\t" "\n"|cat -n|grep ~{pos_col}|cut -f1|tr -d " ")
    #tabix the file
    tabix -s $chrnum -b $posnum -e $posnum ~{sumstats}
    python3 /scripts/filter_hits_regions.py --sumstats ~{sumstats} --regions ~{region} \
    --pheno ~{pheno} --pval_threshold ~{mlogp_threshold} \
    --pos_col ~{pos_col} --chr_col ~{chr_col} --ref_col ~{ref_col} --alt_col ~{alt_col} --mlogp_col ~{mlogp_col}  --chroms ~{sep=" " chroms} --out ./ --log info
  >>>
  
  output {
    Array[File] pheno_chrom_regions = glob("*sig_hits_*")
    File gw_sig_res = outfile
  }
  
  runtime {
    cpu: "1"
    docker: "${docker}"
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: "1"
  }
}


task merge_regions {

  input {
    Array[File] hits
    String docker
  }

  String outfile = "regions.txt"

  command <<<
    while read f; do cat $f >> ~{outfile}; done <~{write_lines(hits)}
  >>>

  output {
    File regions = outfile
  }
  
  runtime {
    cpu: "1"
    docker: "${docker}"
    memory: "2 GB"
    disks: "local-disk 2 HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: "1"
  }
}
