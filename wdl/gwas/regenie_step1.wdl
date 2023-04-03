task step1 {

    Array[String] phenolist
    Boolean is_binary
    File grm_bed
    File grm_bim = sub(grm_bed, "\\.bed$", ".bim")
    File grm_fam = sub(grm_bed, "\\.bed$", ".fam")
    String prefix = basename(grm_bed, ".bed")
    File cov_pheno
    String covariates
    Int bsize
    String options

    String docker
    Boolean auto_remove_sex_covar
    String sex_col_name

    Int covariate_inclusion_threshold

    command <<<

  #      set -euxo pipefail

        n_cpu=`grep -c ^processor /proc/cpuinfo`
        is_single_sex=$(zcat ${cov_pheno} | awk -v sexcol=${sex_col_name} -v phenocols=${sep="," phenolist} \
        ' BEGIN{ is_single=1}
          NR==1{
            for(i=1;i<=NF;i++) {h[$i]=i;};
            split(phenocols,ps, ",");
            prev="";
            is_single=1;
            for(pi in ps) {
              if(!(ps[pi] in h)) {
                 print "Given phenotype " ps[pi] " does not exist in phenofile" > "/dev/stderr"; exit 1;
              }
            }
            if(!(sexcol in h)) {
              print "Given sexcolumn:"sexcol" does not exist in phenofile" > "/dev/stderr"; exit 1;
            }
          }
          NR>1{
            for(pi in ps) {
             if ($h[ps[pi]]!="NA") {
               if(prev!=""&&prev!=$h[sexcol]) {is_single=0; exit;};
               prev=$h[sexcol];
             }
            }
         }
         END { print "is single ", is_single > "/dev/stderr"; printf is_single }'
        )

      if [[ $is_single_sex -eq 1 ]]
      then
        echo "true" > is_single_sex
      else
        echo "false" > is_single_sex
      fi

       if [[ "${auto_remove_sex_covar}" == "true" && "$is_single_sex" == "1" ]];
       then
          covars=$(echo ${covariates} | sed -e 's/${sex_col_name}//' | sed 's/^,//' | sed -e 's/,$//' | sed 's/,,/,/g')
       else
          covars=${covariates}
       fi
        zcat ${cov_pheno}|awk '
        BEGIN {
        FS = "\t"
        _items = split("${sep="," phenolist}",phenos,",")
        }
        NR == 1 {
        for(i = 1; i <= NF; i++) {
        h[$i] = i
        }
        for (phkey in phenos){
        exists=phenos[phkey] in h
        if (!exists) {
        print "Phenotype:"phenos[phkey]" not found in the given phenotype file." > "/dev/stderr"
        err = 1
        exit 1
        }
        }

        }
        NR > 1 {
        for (pk in phenos){
        if ($h[phenos[pk]] != "NA"){
        vals[phenos[pk]$h[phenos[pk]]] +=1
        }
        }
        }
        END{
        if (!err) {
        for (pk in phenos){
        printf "{\"phenotype\":\"%s\",\"cases\":%d,\"controls\":%d}",phenos[pk],vals[phenos[pk]"1"],vals[phenos[pk]"0"] > phenos[pk]"_cases_controls.json"
        }
        }
        }
        '

        #filter out covariates with too few observations
        COVARFILE=${cov_pheno}
        PHENO="${sep=',' phenolist}"
        THRESHOLD=${covariate_inclusion_threshold}
        # Filter binary covariates that don't have enough covariate values in them
        # Inputs: covariate file, comma-separated phenolist, comma-separated covariate list, threshold for excluding a covariate
        # If a covariate is quantitative (=having values different from 0,1,NA), it is masked and will be passed through.
        # If a binary covariate has value NA, it will not be counted towards 0 or 1 for that covariate.
        # If a covariate does not seem to exist (e.g. PC{1:10}), it will be passed through.
        # If any of the phenotypes is not NA on that row, that row will be counted. This is in line with the step1 mean-imputation for multiple phenotypes.
        zcat $COVARFILE |awk -v covariates=$covars  -v phenos=$PHENO -v th=$THRESHOLD > new_covars  '
        BEGIN{FS="\t"}
        NR == 1 {
            covlen = split(covariates,covars,",")
            phlen = split(phenos,phenoarr,",")
            for (i=1; i<=NF; i++){
                h[$i] = i
                mask[$i] = 0
            }
        }
        NR > 1 {
            #if any of the phenotypes is not NA, then take the row into account
            process_row=0
            for (ph in phenoarr){
                if ($h[phenoarr[ph]] != "NA"){
                    process_row = 1
                }
            }
            if (process_row == 1){
                for (co in covars){
                    if($h[covars[co]] == 0) {
                        zerovals[covars[co]] +=1
                    }
                    else if($h[covars[co]] == 1) {
                        onevals[covars[co]] +=1
                    }
                    else if($h[covars[co]] == "NA"){
                        #no-op
                        na=0;
                    }
                    else {
                        #mask this covariate to be included, no matter the counts
                        #includes both covariate not found in header and quantitative covars
                        mask[covars[co]] =1
                    }
                }
            }

        }
        END{
            SEP=""
            for (co in covars){
                if( ( zerovals[covars[co]] > th && onevals[covars[co]] > th ) || mask[covars[co]] == 1 ){
                    printf("%s%s" ,SEP,covars[co])
                    SEP=","
                }
                printf "Covariate %s zero count: %d one count: %d mask: %d\n",covars[co],zerovals[covars[co]],onevals[covars[co]],mask[covars[co]] >> "/dev/stderr";
            }
        }'

        n_cpu=`grep -c ^processor /proc/cpuinfo`
        NEWCOVARS=$(cat new_covars)
        # fid needs to be the same as iid in fam
        awk '{$1=$2} 1' ${grm_fam} > tempfam && mv tempfam ${grm_fam}

        regenie \
        --step 1 \
        ${if is_binary then "--bt" else ""} \
        --bed ${sub(grm_bed, "\\.bed$", "")} \
        --covarFile ${cov_pheno} \
        --covarColList $NEWCOVARS \
        --phenoFile ${cov_pheno} \
        --phenoColList ${sep="," phenolist} \
        --bsize ${bsize} \
        --lowmem \
        --lowmem-prefix tmp_rg \
        --gz \
        --threads $n_cpu \
        --out ${prefix} \
        ${if is_binary then "--write-null-firth" else ""} \
        ${options}

        # rename loco files with phenotype names and update pred.list accordingly giving it a unique name
        awk '{orig=$2; sub(/_[0-9]+.loco.gz/, "."$1".loco.gz", $2); print "mv "orig" "$2} ' ${prefix}_pred.list | bash
        phenohash=`echo ${sep="," phenolist} | md5sum | awk '{print $1}'`
        awk '{sub(/_[0-9]+.loco.gz/, "."$1".loco.gz", $2)} 1' ${prefix}_pred.list > ${prefix}.$phenohash.pred.list
        loco_n=$(wc -l ${prefix}.$phenohash.pred.list|cut -d " " -f 1)

        if [[ "${is_binary}" == "true" ]]
        then
            # rename firth files with phenotype names and update firth.list accordingly giving it a unique name
            awk '{orig=$2; sub(/_[0-9]+.firth.gz/, "."$1".firth.gz", $2); print "mv "orig" "$2} ' ${prefix}_firth.list | bash
            phenohash=`echo ${sep="," phenolist} | md5sum | awk '{print $1}'`
            awk '{sub(/_[0-9]+.firth.gz/, "."$1".firth.gz", $2)} 1' ${prefix}_firth.list > ${prefix}.$phenohash.firth.list

            #check that firth nulls were created
            firth_n=$(wc -l ${prefix}.$phenohash.firth.list|cut -d " " -f 1)
            #check if there is a firth approx per every null
            if [ "$loco_n" -ne "$firth_n" ]; then
                echo "fitting firth null approximations FAILED. This job will abort."
                exit 1
            fi
        else # touch files to have quant phenos not fail output globbing
            touch ${prefix}.$phenohash.firth.list
            touch get_globbed.firth.gz
        fi
    >>>

    output {
        File log = prefix + ".log"
        Array[File] loco = glob("*.loco.gz")
        File pred = glob("*.pred.list")[0]
        Array[File] nulls = glob("*.firth.gz")
        File firth_list = glob("*.firth.list")[0]
        String covars_used = read_string("new_covars")
        File covariatelist = "new_covars"
        Boolean  is_single_sex = read_boolean("is_single_sex")
        Array[File] case_control_counts = glob("*_cases_controls.json")
    }

    runtime {

        docker: "${docker}"
        cpu: if length(phenolist) == 1 then 1 else if length(phenolist) <=10 then 2 else 4
        memory: if length(phenolist) == 1 then "8 GB" else "10 GB"
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}

workflow regenie_step1 {

    Array[String] phenolist
    Boolean is_binary
    String cov_pheno
    String covariates
    String docker
    Boolean auto_remove_sex_covar
    String sex_col_name

    call step1 {
        input: phenolist=phenolist, is_binary=is_binary, cov_pheno=cov_pheno, covariates=covariates,
        auto_remove_sex_covar=auto_remove_sex_covar, docker=docker, sex_col_name=sex_col_name
    }

    output {
        File pred = step1.pred
        Array[File] loco = step1.loco
        File firth_list = step1.firth_list
        Array[File] nulls = step1.nulls
        String covars_used = step1.covars_used
        Boolean is_single_sex = step1.is_single_sex
        Array[File] case_control_counts = step1.case_control_counts
    }
}
