task step2 {

    Array[String] phenolist
    File cov_pheno
    String covariates
    String test
    # regenie accepts either no test specified (i.e. additive), or '--test dominant' or '--test recessive'. '--test additive' is an error.
    String test_cmd = if test == "additive" then "" else "--test "+ test
    Boolean is_binary
    File bgen
    File bgi = bgen + ".bgi"
    File sample = bgen + ".sample"
    File pred
    File firth_list
    String prefix = sub(basename(pred), "_pred.list", "") + "." + basename(bgen)
    Array[File] loco
    Array[File] nulls
    Int bsize
    String options

    String docker

    command <<<

        n_cpu=`grep -c ^processor /proc/cpuinfo`

        # move loco files to /cromwell_root as pred file paths point there
        for file in ${sep=" " loco}; do
            mv $file .
        done

        # move null files to /cromwell_root as firth_list file paths point there
        for file in ${sep=" " nulls}; do
            mv $file .
        done

        regenie \
        --step 2 \
        ${test_cmd} \
        ${if is_binary then "--bt --af-cc" else ""} \
        --bgen ${bgen} \
        --ref-first \
        --sample ${sample} \
        --covarFile ${cov_pheno} \
        --covarColList ${covariates} \
        --phenoFile ${cov_pheno} \
        --phenoColList ${sep="," phenolist} \
        --pred ${pred} \
        ${if is_binary then "--use-null-firth ${firth_list}" else ""} \
        --bsize ${bsize} \
        --threads $n_cpu \
        --gz \
        --out ${prefix} \
        ${options}

    >>>

    output {
        File log = prefix + ".log"
        Array[File] regenie = glob("*.regenie.gz")
    }

    runtime {
        docker: "${docker}"
        cpu: if length(phenolist) == 1 then 1 else if length(phenolist) <=4 then 2 else if length(phenolist) <= 10 then 4 else if length(phenolist) < 16 then 8 else 16
        memory: if length(phenolist) <= 2 then "4 GB" else "6 GB"
        disks: "local-disk " + (ceil(size(bgen, "G")) + 5) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}

task gather {

    Boolean is_binary
    Array[File] files

    String docker

    command <<<

        set -euxo pipefail

        pheno=`basename ${files[0]} .regenie.gz | awk -F "." '{sub(/[^_]*_/, "", $NF); print $NF}'`
        mkdir regenie munged

        echo -e "`date`\tconcatenating result pieces into regenie/$pheno.regenie.gz, sorting by chr pos just in case"
        cat \
        <(zcat ${files[0]} | head -1) \
        <(for file in ${sep=" " files}; do
            zcat $file | tail -n+2
        done | sort -k1,1g -k2,2g) | bgzip > regenie/$pheno.regenie.gz

        if [[ "${is_binary}" == "true" ]]; then
            echo -e "`date`\tconverting to munged/$pheno.gz to a format used for importing to pheweb, omitting variants with -log10p NA (unsuccessful Firth/SPA)"
            zcat regenie/$pheno.regenie.gz | awk '
            BEGIN {FS=" "; OFS="\t"; split("CHROM GENPOS ALLELE0 ALLELE1 LOG10P BETA SE A1FREQ A1FREQ_CASES A1FREQ_CONTROLS", REQUIRED_FIELDS)}
            NR==1 {for(i=1;i<=NF;i++) h[$i]=i;
                   for(i in REQUIRED_FIELDS) if (!(REQUIRED_FIELDS[i] in h)) {print REQUIRED_FIELDS[i]" expected in regenie header">>"/dev/stderr"; exit 1}
                   print "#chrom","pos","ref","alt","pval","mlogp","beta","sebeta","af_alt","af_alt_cases","af_alt_controls"}
            NR>1 && $h["LOG10P"]!="NA" {print $h["CHROM"],$h["GENPOS"],$h["ALLELE0"],$h["ALLELE1"],10^-$h["LOG10P"],$h["LOG10P"],$h["BETA"],$h["SE"],$h["A1FREQ"],$h["A1FREQ_CASES"],$h["A1FREQ_CONTROLS"]}' \
            | bgzip > munged/$pheno.gz
        else
            echo -e "`date`\tconverting to munged/$pheno.gz to a format used for importing to pheweb"
            zcat regenie/$pheno.regenie.gz | awk '
            BEGIN {FS=" "; OFS="\t"; split("CHROM GENPOS ALLELE0 ALLELE1 LOG10P BETA SE A1FREQ", REQUIRED_FIELDS)}
            NR==1 {for(i=1;i<=NF;i++) h[$i]=i;
                   for(i in REQUIRED_FIELDS) if (!(REQUIRED_FIELDS[i] in h)) {print REQUIRED_FIELDS[i]" expected in regenie header">>"/dev/stderr"; exit 1}
                   print "#chrom","pos","ref","alt","pval","mlogp","beta","sebeta","af_alt"}
            NR>1  {print $h["CHROM"],$h["GENPOS"],$h["ALLELE0"],$h["ALLELE1"],10^-$h["LOG10P"],$h["LOG10P"],$h["BETA"],$h["SE"],$h["A1FREQ"]}' \
            | bgzip > munged/$pheno.gz
        fi

        echo -e "`date`\tprinting only chr, pos, pval to speed up and reduce memory use of qqplot.R"
        zcat munged/$pheno.gz | cut -f1,2,5 > $pheno

        echo -e "`date`\trunning qqplot.R"
        qqplot.R --file $pheno --chrcol "#chrom" --bp_col "pos" --pval_col "pval" --loglog_pval 10 > qq_out 2> qq_err

        echo -e "`date`\ttabixing munged/$pheno.gz"
        tabix -s1 -b2 -e2 munged/$pheno.gz
        echo -e "`date`\tdone"

    >>>

    output {
        File regenie = glob("regenie/*.regenie.gz")[0]
        File pheweb = glob("munged/*.gz")[0]
        File pheweb_tbi = glob("munged/*.gz.tbi")[0]
        File qq_out = "qq_out"
        File qq_err = "qq_err"
        Array[File] pngs = glob("*.png")
        Array[File] quantiles = glob("*qquantiles.txt")
    }

    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "6 GB"
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}

workflow regenie_step2 {

    String docker
    Array[String] phenolist
    Boolean is_binary
    String cov_pheno
    String covariates

    String pred
    Array[String] loco
    String firth_list
    Array[String] nulls

    File bgenlist
    Array[String] bgens = read_lines(bgenlist)

    scatter (bgen in bgens) {
        call step2 {
            input: docker=docker, phenolist=phenolist, is_binary=is_binary, cov_pheno=cov_pheno, covariates=covariates, pred=pred, loco=loco, nulls=nulls,firth_list=firth_list,bgen=bgen
        }
    }

    Array[Array[String]] pheno_results = transpose(step2.regenie)
    scatter (pheno_result in pheno_results) {
        call gather {
            input: is_binary=is_binary, files=pheno_result
        }
    }

    output {
        Array[File] pheweb = gather.pheweb
        Array[File] pheweb_tbi = gather.pheweb_tbi
        Array[File] qq_out = gather.qq_out
        Array[File] qq_err = gather.qq_err
        Array[Array[File]] pngs = gather.pngs
        Array[Array[File]] quantiles = gather.quantiles
    }
}
