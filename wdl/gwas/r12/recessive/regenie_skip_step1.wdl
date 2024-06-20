import "regenie_step2_only_sub.wdl" as sub

task coding_gather {

    Array[Array[File]] files
    Array[File] files_flat = flatten(files)

    String docker

    command <<<

    cat <(head -n1 ${files_flat[0]}) <(awk 'FNR>1' ${sep=" " files_flat}) | bgzip > coding_variants.txt.gz

    >>>

    output {
        File coding_variants = "coding_variants.txt.gz"
    }

    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 20 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}

task recalculate_covars{

    Array[String] phenolist
    Boolean is_binary
    File cov_pheno
    String covariates
    String docker
    Boolean auto_remove_sex_covar
    String sex_col_name
    String locotemplate
    String nulltemplate

    Int covariate_inclusion_threshold

    command <<<
        set -eux

        n_cpu=`grep -c ^processor /proc/cpuinfo`
        is_single_sex=$(zcat ${cov_pheno} | awk -v sexcol=${sex_col_name} -v phenocols=${sep="," phenolist} \
        ' BEGIN{ is_single=1}
          NR==1{
            for(i=1;i<=NF;i++) {h[$i]=i;};
            split(phenocols,ps, ",");
            prev="";
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
        zcat $COVARFILE | awk -v covariates=$covars  -v phenos=$PHENO -v th=$THRESHOLD > new_covars  '
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

        for p in ${sep=" " phenolist};do
            echo  "${nulltemplate}"|sed "s/PHENO/$p/" >> nulls
            echo  "${locotemplate}"|sed "s/PHENO/$p/" >> locos
        done
    >>>

    output{
        String covars = read_string("new_covars")
        Array[String] locos = read_lines("locos")
        Array[String] nulls = read_lines("nulls")
    }
    
    runtime{
        docker: "${docker}"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 50 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}

workflow regenie {

    File phenolist
    String covariates
    Array[Array[String]] phenos = read_tsv(phenolist)
    Boolean is_binary
    String cov_pheno
    String nullsfiletemplate
    String locofiletemplate

    scatter (i in range(length(phenos))) {

        Array[String] pheno = phenos[i]
        call recalculate_covars{
            input:
            phenolist=pheno,
            is_binary=is_binary,
            cov_pheno=cov_pheno,
            covariates =covariates,
            locotemplate = locofiletemplate,
            nulltemplate = nullsfiletemplate
        }
        
        Array[String] loco = recalculate_covars.locos
        Array[String] nulls = recalculate_covars.nulls
        
        call sub.regenie_step2 as sub_step2 {
            input: phenos=pheno,
            is_binary=is_binary,
            cov_pheno=cov_pheno,
            covariates=recalculate_covars.covars,
            loco=recalculate_covars.locos,
            nulls=recalculate_covars.nulls
        }
    }

    call coding_gather {
        input: files=sub_step2.coding
    }
}