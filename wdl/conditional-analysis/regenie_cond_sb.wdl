version 1.0

workflow regenie_cond_sb {

  input {

    # FILES & RELEVANT DATA
    String pheno
    String sumstats_root
    File phenofile

    # RELEASE SPECIFIC STUFF
    String docker
    String prefix

    # PARAMS FOR EXTRACTING HITS
    Float locus_mlogp_threshold
    Float conditioning_mlogp_threshold

    Array[String] chroms
    #COLUMN NAMES ET AL
    String chr_col
    String pos_col
    String ref_col
    String alt_col
    String p_col
    String mlogp_col
    String beta_col
    String se_col
    String freq_col
  }
  
  File sumstats = sub(sumstats_root,"PHENO",pheno)
  

  # call that keeps only covariates with minimum count, else regene might fail
  call filter_covariates {input: docker=docker,phenofile=phenofile,pheno=pheno}  

  # this call returns a bed file with all the relevant candidate regions
  call finemap_regions {
    input: pheno=pheno,phenofile=phenofile,sumstats=sumstats, chr_col=chr_col,pos_col=pos_col,ref_col=ref_col,alt_col=alt_col,p_col=p_col,beta_col=beta_col,se_col=se_col,freq_col=freq_col
  }
  #  call that extracts hits > mlogp from the finemap regions generated above
  call extract_cond_regions {
    input: pheno=pheno, mlogp_threshold = locus_mlogp_threshold, docker=docker,mlogp_col=mlogp_col,chr_col=chr_col,pos_col=pos_col,ref_col=ref_col,alt_col=alt_col,chroms=chroms,sumstats=sumstats,region=finemap_regions.bed
  }

  # read in ouotputs of previous steps (covariates and regions) to allow looping
  Map[String,String] cov_map = read_map(filter_covariates.cov_pheno_map)
  Array[Array[String]] all_regions = read_tsv(extract_cond_regions.gw_sig_res)

  # now that we have all valid regions, let's scatter over them so that each VM gets a region
  scatter (region in all_regions) {
    String chrom = region[1]
    String locus_region = region[2] + " " + region[3]
    call regenie_conditional {
      input: docker = docker, prefix=prefix,locus_region=locus_region,pheno=pheno,chrom=chrom,covariates = cov_map[pheno],mlogp_col = mlogp_col,chr_col=chr_col,pos_col = pos_col,ref_col=ref_col,alt_col=alt_col,beta_col=beta_col,se_col = se_col,pval_threshold=conditioning_mlogp_threshold,sumstats_root=sumstats_root,phenofile=phenofile
    }
  }

  # flatten array of arrays into single list of results and merge them into individual chain and log files
  Array[File] results = flatten(regenie_conditional.conditional_chains)
  Array[File] logs = flatten(regenie_conditional.logs)
  call merge_results{input:pheno=pheno,docker=docker,prefix=prefix,result_list = results,logs=logs}
  
}


task merge_results {
  
  input {
    String pheno
    Array[File] result_list
    Array[File] logs
    String docker
    String prefix
  }
  
  String out_file = prefix + "_" + pheno + ".independent.snps.txt"
  String log_file = sub(out_file,".txt",".log")
  
  command <<<
  cat ~{write_lines(result_list)} > results.txt
  cat ~{write_lines(logs)} > logs.txt

  # keep header from one of the files
  head -n1 results.txt | xargs head -n 1 > ~{out_file}
  
  while read f; do cat $f | sed -E 1d >> ~{out_file} ; done < results.txt
  while read f; do cat $f                 >> ~{log_file} ; done < logs.txt
  >>>
  
  output {
    File pheno_independent_snps = out_file
    File pheno_logs = log_file
  }
  
  runtime {
    cpu: "4"
    docker: "${docker}"
    memory: "4 GB"
    disks: "local-disk 10 HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: "1"   
  }
  
}
  
  
task regenie_conditional {
  
  input {
    # GENERAL PARAMS
    String docker
    String prefix
    # hit info
    String locus_region
    String pheno
    String chrom
    # files to localize 
    File phenofile
    String bgen_root
    String null_root
    String sumstats_root
    # column names and stuff
    String chr_col
    String pos_col
    String ref_col
    String alt_col
    String mlogp_col
    String beta_col
    String se_col
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
  Int disk_size = ceil(size(bgen,'GB')) + ceil(size(sumstats,'GB')) + ceil(size(null,'GB')) + ceil(size(phenofile,'GB')) + 1


  command <<<
    
    echo ~{pheno} ~{chrom} ~{cpus} 
    tabix -h ~{sumstats}  ~{chrom} > region_sumstats.txt

    python3 /scripts/regenie_conditional.py \
    --out ./~{prefix}  --bgen ~{bgen}  --null-file ~{null}  --sumstats region_sumstats.txt \
    --pheno-file ~{phenofile} --pheno ~{pheno} \
    --locus-region ~{locus_region}  --pval-threshold ~{pval_threshold} --max-steps ~{max_steps} \
    --chr-col "~{chr_col}" --pos-col "~{pos_col}" --ref-col "~{ref_col}" --alt-col "~{alt_col}" --mlogp-col "~{mlogp_col}" --beta-col "~{beta_col}" --sebeta-col "~{se_col}" \
    --covariates ~{covariates} ~{if defined(regenie_params) then " --regenie-params " + regenie_params else ""} --log info

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

  File index = sumstats + ".tbi"
  Int disk_size = ceil(size(sumstats,'GB')) + ceil(size(region,'GB')) + 1
  
  String outfile= pheno + "_sig_hits.txt" 
  
  command <<<

    python3 /scripts/filter_hits_regions.py --sumstats ~{sumstats} --regions ~{region} \
    --pheno ~{pheno} --pval_threshold ~{mlogp_threshold} \
    --pos_col "~{pos_col}" --chr_col "~{chr_col}" --ref_col "~{ref_col}" --alt_col "~{alt_col}" --mlogp_col "~{mlogp_col}"  --chroms ~{sep=" " chroms} --out ./ --log info
  >>>
  
  output {
    Array[File] pheno_chrom_regions = glob("*sig_hits_*")
    File gw_sig_res = outfile
  }
  
  runtime {
    cpu: "1"
    docker: "${docker}"
    memory: "2 GB"
    disks: "local-disk ${disk_size} HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: "1"
  }
}




task filter_covariates {

  input {
    File phenofile
    Array[String] covariates
    String pheno
    String docker
    Int threshold_cov_count
    }

    String outfile = "./pheno_cov_map_" + threshold_cov_count + ".txt"
    Int disk_size = ceil(size(phenofile,'GB')) + 2 * 2
    
    command <<<

      set -euxo pipefail
      echo ~{pheno} > pheno_list.txt
      python3 <<CODE
      
      import pandas as pd
      import numpy as np
      
      #read in phenos as list of phenos regardless
      tot_phenos = []
      phenos_groups = []
      with open('pheno_list.txt') as i:
          for line in i:
              phenos = line.strip().split()
              phenos_groups.append(phenos)
              tot_phenos += phenos    

      #read in phenos mapping all valid entries to 1 and NAs to 0
      pheno_df= pd.read_csv('~{phenofile}',sep='\t',usecols=tot_phenos).notna().astype(int)
      print(pheno_df)
      # read in covariates getting absolute values (handles PCs)
      covariates= '~{sep="," covariates}'.split(',')
      cov_df= pd.read_csv('~{phenofile}',sep='\t',usecols=covariates).abs()
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



task finemap_regions {

  input {
    # THINGS THAT ARE INHERITED FROM GLOBAL PARAMS
    String pheno
    File phenofile
    File sumstats
    String chr_col
    String pos_col
    String ref_col
    String alt_col
    String p_col
    String freq_col
    String beta_col
    String se_col

    #RUNTIME
    String docker
    String zones
    Int cpu
    Int mem

    #OTHER STUFF
    Boolean scale_se_by_pval
    Boolean x_chromosome
    Boolean set_variant_id
    String rsid_col
    String delimiter
    Int window
    Int max_region_width
    Float window_shrink_ratio
    # can be helpful if adding finemapping with relaxed threshold after more stringent has already ben run.
    # does not include regions with lead snp < this
    Float p_threshold
    Float? minimum_pval
    String? set_variant_id_map_chr
  }
  
  
    
  command <<<
  
  catcmd="cat"
  if [[ ~{phenofile} == *.gz ]] || [[ ~{phenofile} == *.bgz ]]
  then
        catcmd="zcat"
  fi
   
  echo "Reading phenotype file with $catcmd"
  $catcmd ~{phenofile} | awk -v ph=~{pheno} '
        BEGIN {
            FS = "\t"
        }
        NR == 1 {
            for(i = 1; i <= NF; i++) {
                h[$i] = i
            }
            exists=ph in h
            if (!exists) {
                print "Phenotype:"ph" not found in the given phenotype file." > "/dev/stderr"
                err = 1
                exit 1
            }
        }
        NR > 1 && $h[ph] != "NA" {
            vals[$h[ph]] += 1
            print $1 > ph".incl"
            if ($h[ph] != 0 && $h[ph] != 1 && !err) {
                print "Phenotype:"ph" seems a quantitative trait. Setting var_y = 1 and prior_std = 0.05." > "/dev/stderr"
                print 1.0 > "var_y.txt"
                print 0.05 > "prior_std.txt"
                err = 1
            }
        }
        END {
            if (!err) {
                phi = vals["1"] / (vals["1"]+vals["0"])
                var_y = phi * (1-phi)
                std = 0.05 * sqrt(phi*(1-phi))
                print var_y > "var_y.txt"
                print std > "prior_std.txt"
            }
        }'
  
  if [[ $? -ne 0 ]]
  then
      echo "Error occurred while getting case control counts for ~{pheno}"
      exit 1
  fi
  
  wc -l ~{pheno}.incl | cut -f1 -d' ' > n_samples.txt
  
  make_finemap_inputs.py \
      --sumstats ~{sumstats} \
      --rsid-col "~{rsid_col}" \
      --chromosome-col "~{chr_col}" \
      --position-col "~{pos_col}" \
      --allele1-col "~{ref_col}" \
      --allele2-col "~{alt_col}" \
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
      
    File bed = pheno + ".bed"
    File log = pheno + ".log"
    
  }
    
  runtime {
    docker: "${docker}"
    cpu: "${cpu}"
    memory: "${mem} GB"
    disks: "local-disk 20 HDD"
    zones: "${zones}"
    preemptible: 2
    noAddress: true
  }
}
