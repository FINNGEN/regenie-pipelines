import "regenie_sub.wdl" as sub

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
        ${options}

    >>>

    output {
        File log = prefix + ".log"
        Array[File] loco = glob("*.loco.gz")
        File pred = prefix + "_pred.list"
    }

    runtime {

        docker: "${docker}"
        cpu: if length(phenolist) == 1 then 1 else if length(phenolist) <=10 then 2 else 4
        memory: if length(phenolist) == 1 then "6 GB" else "8 GB"
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    }
}

workflow regenie {

    File phenolist
    Array[Array[String]] pheno_chunks = read_tsv(phenolist)
    Boolean is_binary
    String cov_pheno
    String covariates

    scatter (pheno_chunk in pheno_chunks) {
        call step1 {
            input: phenolist=pheno_chunk, is_binary=is_binary, cov_pheno=cov_pheno, covariates=covariates
        }
        call sub.regenie_step2 as sub_step2 {
            input: phenolist=pheno_chunk, is_binary=is_binary, cov_pheno=cov_pheno, covariates=covariates, pred=step1.pred, loco=step1.loco
        }
    }
}
