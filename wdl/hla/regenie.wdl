import "regenie_step1.wdl" as step1
import "regenie_sub.wdl" as sub

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

}
