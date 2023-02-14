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

    Boolean auto_remove_sex_covar
    String sex_col_name

    scatter (pheno_chunk in pheno_chunks) {
        call step1.regenie_step1 as sub_step1 {
            input: phenolist=pheno_chunk, is_binary=is_binary, cov_pheno=cov_pheno, covariates=covariates,
            auto_remove_sex_covar=auto_remove_sex_covar,sex_col_name=sex_col_name
        }
        call sub.regenie_step2 as sub_step2 {
            input: phenolist=pheno_chunk, is_binary=is_binary, cov_pheno=cov_pheno, covariates=sub_step1.covars_used,
            pred=sub_step1.pred, loco=sub_step1.loco, nulls=sub_step1.nulls, firth_list=sub_step1.firth_list,
            sex_col_name=sex_col_name, is_single_sex=sub_step1.is_single_sex
        }
    }

    call coding_gather {
        input: files=sub_step2.coding
    }

output
                    {
                       
                        File coding_var=coding_gather.coding_variants
                        Array[File] pheweb_sumst=flatten(sub_step2.pheweb)
                        Array[File] pheweb_sumst_tbi=flatten(sub_step2.pheweb_tbi)
                        Array[File] plots=flatten(flatten(sub_step2.pngs))
                        Array[File] quant=flatten(flatten(sub_step2.quantiles))
                        Array[File] sig_stat=flatten(sub_step2.summary)
                        Array[File] casecontNjson=flatten(sub_step1.case_control_counts) 
                    }
}
