import "base_sub.wdl" as sub

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
    String cpuplatform

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

        NEWCOVARS=$(cat new_covars)

        for p in ${sep=" " phenolist};do
            echo  "${nulltemplate}"|sed "s/PHENO/$p/" >> nulls
            echo  "${locotemplate}"|sed "s/PHENO/$p/" >> locos
        done

        ## calculate sex-specific samples
        zcat ${cov_pheno} | awk 'NR==1{ for(i=1;i<=NF;i++) {h[$i]=i};
                                            if(!("${sex_col_name}" in h)) {
                                            print "Given sex column not found in phenotype file" > "/dev/stderr";
                                            exit 1;
                                            }
                                    }
                                NR>1&&$h["${sex_col_name}"]==0{ print $1,$1}' > males

        zcat ${cov_pheno} | awk 'NR==1{  for(i=1;i<=NF;i++) {h[$i]=i };
                                        if(!("${sex_col_name}" in h))
                                        {
                                        print "Given sex column not found in phenotype file" > "/dev/stderr";
                                        exit 1;
                                        }
                                    }
                                NR>1&&$h["${sex_col_name}"]==1{ print $1,$1}' > females

        sex_covars=$(echo $NEWCOVARS | sed -e 's/${sex_col_name}//' | sed 's/^,//' | sed -e 's/,$//' | sed 's/,,/,/g')
        echo $sex_covars > sex_edited_covars

        if [[ $sex_covars == $NEWCOVARS ]];
        then
        echo "Warning! No sex covariate detected in used covariates."
        fi
        for p in ${sep=" " phenolist}; do
            if [[ "${is_binary}" == "true" ]];
            then
                N_cases_females=$(awk -v pheno=$p 'FNR==NR { females[$1]; next } FNR==1{ for(i=1;i<=NF;i++) {h[$i]=i} }; NR>1 && $h[pheno]==1 && $1 in females {print $1}' females <(zcat ${cov_pheno}) | wc -l)
                N_cases_males=$(awk -v pheno=$p 'FNR==NR { males[$1]; next } FNR==1{ for(i=1;i<=NF;i++) {h[$i]=i} }; NR>1 && $h[pheno]==1 && $1 in males {print $1}' males <(zcat ${cov_pheno}) | wc -l)
            else 
                N_cases_females=$(awk -v pheno=$p 'FNR==NR { females[$1]; next } FNR==1{ for(i=1;i<=NF;i++) {h[$i]=i} }; NR>1 && $h[pheno]!="NA" && $1 in females {print $1}' females <(zcat ${cov_pheno}) | wc -l)
                N_cases_males=$(awk -v pheno=$p 'FNR==NR { males[$1]; next } FNR==1{ for(i=1;i<=NF;i++) {h[$i]=i} }; NR>1 && $h[pheno]!="NA" && $1 in males {print $1}' males <(zcat ${cov_pheno}) | wc -l)
            fi 
            if [[ $N_cases_females -lt 10 ]] || [[ $N_cases_males -lt 10 ]];
            then
                echo "Less than 10 cases in a sex. Skipping testing of sex specific effects for endpoint $p."
                continue
            fi
            echo $p >> endpoints_to_process
        done
        # test if there are any phenos
        if test -f endpoints_to_process; then
            echo "true" > analyse_endpoints
        else
            touch endpoints_to_process
            echo "false" > analyse_endpoints
        fi 
    >>>

    output{
        String covars = read_string("sex_edited_covars")
        Array[String] locos = read_lines("locos")
        Array[String] nulls = read_lines("nulls")
        Boolean do_analysis = read_boolean("analyse_endpoints")
        Array[String] phenotypes_to_process = if do_analysis then read_lines("endpoints_to_process") else []
        File males = "males"
        File females = "females"
        
    }
    
    runtime{
        docker: "${docker}"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 50 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
        #cpuPlatform:cpuplatform
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
    String cpuplatform

    scatter (i in range(length(phenos))) {

        Array[String] pheno = phenos[i]
        call recalculate_covars{
            input:
            phenolist=pheno,
            is_binary=is_binary,
            cov_pheno=cov_pheno,
            covariates =covariates,
            locotemplate = locofiletemplate,
            nulltemplate = nullsfiletemplate,
            cpuplatform=cpuplatform
        }
        
        Array[String] loco = recalculate_covars.locos
        Array[String] nulls = recalculate_covars.nulls
        if(recalculate_covars.do_analysis){
            call sub.regenie_step2 as sub_step2 {
                input: phenos=recalculate_covars.phenotypes_to_process,
                is_binary=is_binary,
                cov_pheno=cov_pheno,
                covariates=recalculate_covars.covars,
                loco=recalculate_covars.locos,
                nulls=recalculate_covars.nulls,
                females = recalculate_covars.females,
                males = recalculate_covars.males,
                cpuplatform=cpuplatform
            }
        }
    }

    output{
        Array[File] pheweb = flatten(select_all(sub_step2.pheweb))
        Array[File] pheweb_tbi = flatten(select_all(sub_step2.pheweb_tbi))
        Array[File] qq_out = flatten(select_all(sub_step2.qq_out))
        Array[File] qq_err = flatten(select_all(sub_step2.qq_err))
        Array[File] pngs = flatten(select_all(sub_step2.pngs))
        Array[File] summaries = flatten(select_all(sub_step2.summaries))
        Array[File] coding = flatten(select_all(sub_step2.coding))
    }
}
