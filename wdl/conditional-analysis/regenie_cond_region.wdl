version 1.0

workflow conditional_analysis {

  input {
    File cond_regions
    String pheno
    String prefix
    String mlogp_col
    String chr_col
    String pos_col
    String ref_col
    String alt_col
    Float conditioning_mlogp_threshold
    Array[String] covariates
    File pheno_file
    String sumstats_root
  }

  String docker = "eu.gcr.io/finngen-refinery-dev/regenie:3.3_cond_region"
  call filter_covariates {input: docker=docker,pheno_file=pheno_file,pheno_list = write_lines([pheno]),covariates=covariates}  
  Map[String,String] cov_map = read_map(filter_covariates.cov_pheno_map)

  Array[Array[String]] locus_data = read_tsv(cond_regions)
  scatter (entry in locus_data) {
    String chrom = entry[0]
    String locus_range = entry[0] + ":" + entry[1]
    String locus = entry[2]
    call regenie_conditional {
      input:
      docker = docker,
      prefix=prefix,
      locus=locus,
      region=locus_range,
      pheno=pheno,chrom=chrom,
      covariates = cov_map[pheno],mlogp_col = mlogp_col,chr_col=chr_col,pos_col = pos_col,ref_col=ref_col,alt_col=alt_col,pval_threshold=conditioning_mlogp_threshold,sumstats_root=sumstats_root,pheno_file=pheno_file
      }    
   }
}


task regenie_conditional {
  
  input {
    # GENERAL PARAMS
    String? regenie_docker
    String docker
    String prefix
    # hit info
    String locus
    String region
    String pheno
    String chrom
    # files to localize 
    File pheno_file
    String bgen_root
    String null_root
    String sumstats_root
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
  File sumstats = sub(sumstats_root,"PHENO",pheno)
  File sum_tabix = sumstats + ".tbi"
  File null =sub(null_root,"PHENO",pheno)
  File bgen = sub(bgen_root,'CHROM',chrom)
  File bgen_sample = bgen + ".sample"
  File bgen_index = bgen + ".bgi"


  # runtime params based on file sizes
  Int disk_size = ceil(size(bgen,'GB')) + ceil(size(sumstats,'GB')) + ceil(size(null,'GB')) + ceil(size(pheno_file,'GB')) + 1
  String final_docker = if defined(regenie_docker) then regenie_docker else docker


  command <<<
    
    echo ~{pheno} ~{chrom} ~{cpus} 
    tabix -h ~{sumstats}  ~{region} > region_sumstats.txt

    python3 /scripts/regenie_conditional.py \
    --out ./~{prefix}  --bgen ~{bgen}  --null-file ~{null}  --sumstats region_sumstats.txt \
    --pheno-file ~{pheno_file} --pheno ~{pheno} \
    --locus-region ~{locus} ~{region}  --pval-threshold ~{pval_threshold} --max-steps ~{max_steps} \
    --chr-col ~{chr_col} --pos-col ~{pos_col} --ref-col ~{ref_col} --alt-col ~{alt_col} --mlogp-col ~{mlogp_col} --beta-col ~{beta} --sebeta-col ~{sebeta} \
    --covariates ~{covariates} ~{if defined(regenie_params) then " --regenie-params " + regenie_params else ""} --log info

  >>>
  output {
    Array[File] conditional_chains = glob("./${prefix}*.snps")
    Array[File] logs = glob("./${prefix}*.log")
    Array[File] regenie_output = glob("./${prefix}*.conditional")    
  }
  
  runtime {
    cpu: "~{cpus}"
    docker: "${final_docker}"
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

