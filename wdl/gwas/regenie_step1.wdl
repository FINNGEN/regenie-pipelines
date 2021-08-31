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

    command <<<

        set -euxo pipefail

        n_cpu=`grep -c ^processor /proc/cpuinfo`

        # fid needs to be the same as iid in fam
        awk '{$1=$2} 1' ${grm_fam} > tempfam && mv tempfam ${grm_fam}

        regenie \
        --step 1 \
        ${if is_binary then "--bt" else ""} \
        --bed ${sub(grm_bed, "\\.bed$", "")} \
        --covarFile ${cov_pheno} \
        --covarColList ${covariates} \
        --phenoFile ${cov_pheno} \
        --phenoColList ${sep="," phenolist} \
        --bsize ${bsize} \
        --lowmem \
        --lowmem-prefix tmp_rg \
        --gz \
        --threads $n_cpu \
        --out ${prefix} \
        --write-null-firth
        ${options}

        # rename loco files with phenotype names and update pred.list accordingly giving it a unique name
        awk '{orig=$2; sub(/_[0-9]+.loco.gz/, "."$1".loco.gz", $2); print "mv "orig" "$2} ' ${prefix}_pred.list | bash
        phenohash=`echo ${sep="," phenolist} | md5sum | awk '{print $1}'`
        awk '{sub(/_[0-9]+.loco.gz/, "."$1".loco.gz", $2)} 1' ${prefix}_pred.list > ${prefix}.$phenohash.pred.list

        # rename firth files with phenotype names and update firth.list accordingly giving it a unique name
        awk '{orig=$2; sub(/_[0-9]+.firth.gz/, "."$1".firth.gz", $2); print "mv "orig" "$2} ' ${prefix}_firth.list | bash
        phenohash=`echo ${sep="," phenolist} | md5sum | awk '{print $1}'`
        awk '{sub(/_[0-9]+.firth.gz/, "."$1".firth.gz", $2)} 1' ${prefix}_firth.list > ${prefix}.$phenohash.firth.list

    >>>

    output {
        File log = prefix + ".log"
        Array[File] loco = glob("*.loco.gz")
        File pred = glob("*.pred.list")[0]
        Array[File] nulls = glob("*.firth.gz")
        File firth_list = glob("*.firth.list")[0]
    }

    runtime {

        docker: "${docker}"
        cpu: if length(phenolist) == 1 then 1 else if length(phenolist) <=10 then 2 else 4
        memory: if length(phenolist) == 1 then "6 GB" else "8 GB"
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
    }
}
