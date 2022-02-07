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
    Int cov_filter_th = 10
    command <<<

        set -euxo pipefail


        #filter degenerate covariates, e.g. sex imputed from sex-specific phenotypes
        COVARFILE=${cov_pheno}
        PHENO="${phenolist sep=','}"
        COVARS="${covariates}"
        THRESHOLD=${cov_filter_th}
        # Filter binary covariates that don't have enough covariate values in them
        # Inputs: covariate file, comma-separated phenolist, comma-separated covariate list, threshold for excluding a covariate
        # If a covariate is quantitative (=having values different from 0,1,NA), it is masked and will be passed through.
        # If a binary covariate has value NA, it will not be counted towards 0 or 1 for that covariate.
        # If a covariate does not seem to exist (e.g. PC{1:10}), it will be passed through.
        # If any of the phenotypes is not NA on that row, that row will be counted. This is in line with the step1 mean-imputation for multiple phenotypes.
        zcat $COVARFILE |awk -v covariates=$COVARS  -v phenos=$PHENO -v th=$THRESHOLD > new_covars  '
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
        else
            touch ${prefix}.$phenohash.firth.list
        fi
    >>>

    output {
        File log = prefix + ".log"
        Array[File] loco = glob("*.loco.gz")
        File pred = glob("*.pred.list")[0]
        Array[File] nulls = glob("*.firth.gz")
        File firth_list = glob("*.firth.list")[0]
        String used_covariates = read_string("new_covars")
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

    call step1 {
        input: phenolist=phenolist, is_binary=is_binary, cov_pheno=cov_pheno, covariates=covariates
    }

    output {
        File pred = step1.pred
        Array[File] loco = step1.loco
        File firth_list = step1.firth_list
        Array[File] nulls = step1.nulls
        String new_covariates = step1.new_covars
    }
}
