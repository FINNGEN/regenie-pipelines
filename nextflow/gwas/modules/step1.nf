// Regenie step 1: fits the null model (whole-genome ridge regression) for a
// chunk of phenotypes against a GRM. Mirrors wdl/gwas/regenie_step1.wdl.

process STEP1 {
    tag "${chunk_id} (${phenolist.size()} phenos)"
    container params.step1_docker

    cpus   { phenolist.size() == 1 ? 1 : (phenolist.size() <= 10 ? 2 : 4) }
    memory { 15.GB }
    disk   200.GB, type: 'pd-standard'

    publishDir "${params.outdir}/step1/${chunk_id}", mode: 'copy'

    input:
    tuple val(chunk_id), val(phenolist), path(grm_bed), path(grm_bim), path(grm_fam), path(cov_pheno)

    output:
    tuple val(chunk_id),
          val(phenolist),
          path("*.${chunk_id}.pred.list"),
          path("*.loco.gz"),
          path("*.${chunk_id}.firth.list"),
          path("*.firth.gz", arity: '0..*'),
          path("new_covars"),
          path("is_single_sex"),
          path("${prefix}.log"), emit: main

    script:
    prefix = grm_bed.getBaseName()
    def is_binary        = params.is_binary
    def binary_flag      = is_binary ? '--bt' : ''
    def write_firth_flag = is_binary ? '--write-null-firth' : ''
    def phenolist_csv    = phenolist.join(',')
    def covariates       = params.covariates
    def bsize            = params.step1_bsize
    def threshold        = params.covariate_inclusion_threshold
    def step1_options    = params.step1_options ?: ''
    def sex_col          = params.sex_col_name
    def auto_remove_sex  = params.auto_remove_sex_covar as boolean
    """
    set -eux

    # 1) Single-sex detection
    IS_SINGLE_SEX=\$(zcat ${cov_pheno} | detect_single_sex.awk -v sexcol=${sex_col} -v phenocols=${phenolist_csv})
    if [[ \$IS_SINGLE_SEX -eq 1 ]]; then
        echo "true"  > is_single_sex
    else
        echo "false" > is_single_sex
    fi

    # 2) Optionally strip the sex column from the covariate list
    if [[ "${auto_remove_sex}" == "true" && "\$IS_SINGLE_SEX" == "1" ]]; then
        covars=\$(echo "${covariates}" | sed -e 's/${sex_col}//' | sed 's/^,//' | sed -e 's/,\$//' | sed 's/,,/,/g')
    else
        covars="${covariates}"
    fi

    # 3) Filter covariates with too few observations
    zcat ${cov_pheno} \\
        | filter_covariates.awk -v covariates="\$covars" -v phenos=${phenolist_csv} -v th=${threshold} \\
        > new_covars
    NEWCOVARS=\$(cat new_covars)

    # 4) Fix FID in fam file (regenie expects FID == IID)
    fix_fam_fid.awk ${grm_fam} > tempfam && mv tempfam ${grm_fam}

    # 5) Run regenie step 1
    n_cpu=\$(grep -c ^processor /proc/cpuinfo)
    regenie \\
        --step 1 \\
        ${binary_flag} \\
        --bed ${prefix} \\
        --covarFile ${cov_pheno} \\
        --covarColList \$NEWCOVARS \\
        --phenoFile ${cov_pheno} \\
        --phenoColList ${phenolist_csv} \\
        --bsize ${bsize} \\
        --lowmem \\
        --lowmem-prefix tmp_rg \\
        --gz \\
        --threads \$n_cpu \\
        --out ${prefix} \\
        ${write_firth_flag} \\
        ${step1_options}

    # 6) Rename loco files from <prefix>_N.loco.gz to <prefix>.<pheno>.loco.gz
    rename_loco_from_predlist.awk ${prefix}_pred.list | bash
    rewrite_predlist.awk ${prefix}_pred.list > ${prefix}.${chunk_id}.pred.list

    loco_n=\$(wc -l < ${prefix}.${chunk_id}.pred.list)
    if [[ \$loco_n -ne ${phenolist.size()} ]]; then
        echo "The model did not converge: expected ${phenolist.size()} loco files, got \$loco_n." >&2
        exit 1
    fi

    # 7) Firth handling
    if [[ "${is_binary}" == "true" ]]; then
        rename_firth_from_firthlist.awk ${prefix}_firth.list | bash
        rewrite_firthlist.awk ${prefix}_firth.list > ${prefix}.${chunk_id}.firth.list
        firth_n=\$(wc -l < ${prefix}.${chunk_id}.firth.list)
        if [[ \$loco_n -ne \$firth_n ]]; then
            echo "fitting firth null approximations FAILED" >&2
            exit 1
        fi
    else
        # Quant phenos: no firth files; create an empty firth.list so outputs match.
        touch ${prefix}.${chunk_id}.firth.list
    fi
    """

    stub:
    prefix = grm_bed.getBaseName()
    def phenolist_csv = phenolist.join(',')
    """
    echo "true" > is_single_sex
    echo "SEX_IMPUTED,AGE" > new_covars
    touch ${prefix}.log
    touch ${prefix}.${chunk_id}.pred.list
    touch ${prefix}.${chunk_id}.firth.list
    ${phenolist.collect { "touch ${prefix}.${it}.loco.gz" }.join('\n    ')}
    ${params.is_binary ? phenolist.collect { "touch ${prefix}.${it}.firth.gz" }.join('\n    ') : ''}
    """
}
