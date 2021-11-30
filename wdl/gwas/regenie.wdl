import "regenie_step1.wdl" as step1
import "regenie_sub.wdl" as sub

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

workflow regenie {

    File phenolist
    Array[Array[String]] pheno_chunks = read_tsv(phenolist)
    Boolean is_binary
    String cov_pheno
    String covariates

    scatter (pheno_chunk in pheno_chunks) {
        call step1.regenie_step1 as sub_step1 {
            input: phenolist=pheno_chunk, is_binary=is_binary, cov_pheno=cov_pheno, covariates=covariates
        }
        call sub.regenie_step2 as sub_step2 {
            input: phenolist=pheno_chunk, is_binary=is_binary, cov_pheno=cov_pheno, covariates=covariates, pred=sub_step1.pred, loco=sub_step1.loco, nulls=sub_step1.nulls, firth_list=sub_step1.firth_list
        }
    }

    call coding_gather {
        input: files=sub_step2.coding
    }
}
