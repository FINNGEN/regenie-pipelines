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

workflow regenie {

    File phenolist
    Array[String] phenos = read_lines(phenolist)
    Boolean is_binary
    String cov_pheno
    String covariates
    String nullsfiletemplate
    String locofiletemplate

    scatter (pheno in phenos) {

        String loco =  sub(locofiletemplate, "PHENO", pheno)
        String nulls =  sub(nullsfiletemplate, "PHENO", pheno)
        
        call sub.regenie_step2 as sub_step2 {
            input: pheno=pheno,
            is_binary=is_binary,
            cov_pheno=cov_pheno,
            covariates=covariates,
            loco=[loco],
            nulls=[nulls]
        }
    }

    call coding_gather {
        input: files=sub_step2.coding
    }
}