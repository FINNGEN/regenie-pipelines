// Regenie step 2: runs association testing per genotype chunk for a pheno
// chunk whose null model was fit in STEP1. One unified process handles both
// --bgen and --pgen inputs — the genotype flags are chosen in Groovy before
// the script block. Mirrors wdl/gwas/regenie_sub.wdl (step2 task) and the
// r13 chip pgen variant.

process STEP2 {
    tag "${chunk_id} x ${geno.getName()}"
    container params.step2_docker

    cpus   { params.step2_resource_profile == 'fixed' ? 1 : (phenolist.size() == 1 ? 1 : (phenolist.size() <= 4 ? 2 : (phenolist.size() <= 10 ? 4 : (phenolist.size() < 16 ? 8 : 16)))) }
    memory { params.step2_resource_profile == 'fixed' ? 6.GB : (phenolist.size() <= 2 ? 4.GB : 6.GB) }
    disk   200.GB, type: 'pd-standard'

    publishDir "${params.outdir}/step2", mode: 'copy', pattern: '*.regenie.gz'
    publishDir "${params.outdir}/step2", mode: 'copy', pattern: '*.sex_spec.gz'

    input:
    tuple val(chunk_id),
          val(phenolist),
          path(pred),
          path(loco_files, stageAs: './*'),
          path(firth_list),
          path(null_files, stageAs: './*'),
          path(cov_pheno),
          val(covars_used),
          val(is_single_sex),
          path(geno),
          path(geno_companion1),
          path(geno_companion2)

    output:
    tuple val(chunk_id),
          val(phenolist),
          val(geno_tag),
          path("${prefix}_*.regenie.gz"),
          path("${prefix}*.sex_spec.gz"),
          path("${prefix}*.log"), emit: main

    script:
    geno_tag = geno.getName()
    prefix   = pred.getName().replaceAll(/_pred\.list$/, '') + '.' + geno_tag
    def is_binary         = params.is_binary
    def binary_flag       = is_binary ? '--bt --af-cc' : ''
    def use_firth_flag    = is_binary ? "--use-null-firth ${firth_list}" : ''
    def test              = params.step2_test
    def test_cmd          = (test == 'additive') ? '' : "--test ${test}"
    def phenolist_csv     = phenolist.join(',')
    def first_pheno       = phenolist[0]
    def bsize             = params.step2_bsize
    def options           = params.step2_options ?: ''
    def sex_col           = params.sex_col_name
    def sex_logp          = params.sex_specific_logpval
    def run_sex_specific  = (params.run_sex_specific as boolean) && !(is_single_sex as boolean)
    def is_binary_flag    = is_binary ? 1 : 0

    def geno_flags
    if (params.input_format == 'pgen') {
        geno_flags = '--pgen ' + geno.getBaseName()
    } else {
        geno_flags = '--bgen ' + geno.getName() + ' --ref-first --sample ' + geno_companion2.getName()
    }

    def do_strip_chr = params.step2_strip_chr as boolean

    """
    set -eux
    n_cpu=\$(grep -c ^processor /proc/cpuinfo)

    # 1) Run regenie step 2 on the whole pheno chunk
    regenie \\
        --step 2 \\
        ${test_cmd} \\
        ${binary_flag} \\
        ${geno_flags} \\
        --covarFile ${cov_pheno} \\
        --covarColList ${covars_used} \\
        --phenoFile ${cov_pheno} \\
        --phenoColList ${phenolist_csv} \\
        --pred ${pred} \\
        ${use_firth_flag} \\
        --bsize ${bsize} \\
        --threads \$n_cpu \\
        --gz \\
        --out ${prefix} \\
        ${options}

    if [[ "${do_strip_chr}" == "true" ]]; then
        # Strip 'chr' prefix from CHROM column (chip pgen variant).
        # The WDL only rewrites the first pheno's file — reproduce that quirk.
        first_out="${prefix}_${first_pheno}.regenie.gz"
        if [[ -s \$first_out ]]; then
            zcat \$first_out | strip_chr_from_chrom.awk | bgzip > \${first_out}.tmp
            mv \${first_out}.tmp \$first_out
        fi
    fi

    # 2) Sex-specific analysis (optional)
    if [[ "${run_sex_specific}" == "true" ]]; then
        zcat ${cov_pheno} | extract_sex_lists.awk -v sexcol=${sex_col} -v sexval=0 > males
        zcat ${cov_pheno} | extract_sex_lists.awk -v sexcol=${sex_col} -v sexval=1 > females

        sex_covars=\$(echo "${covars_used}" | sed -e 's/${sex_col}//' | sed 's/^,//' | sed -e 's/,\$//' | sed 's/,,/,/g')
        if [[ "\$sex_covars" == "${covars_used}" ]]; then
            echo "Warning! No sex covariate detected in used covariates."
        fi

        for p in ${phenolist.join(' ')}; do
            N_cases_females=\$(count_sex_cases.awk -v pheno=\$p -v is_binary=${is_binary_flag} females <(zcat ${cov_pheno}) | wc -l)
            N_cases_males=\$(count_sex_cases.awk -v pheno=\$p -v is_binary=${is_binary_flag} males <(zcat ${cov_pheno}) | wc -l)

            echo "Female cases: \$N_cases_females"
            echo "Male cases:   \$N_cases_males"

            if [[ \$N_cases_females -lt 10 ]] || [[ \$N_cases_males -lt 10 ]]; then
                echo "Less than 10 cases in a sex. Skipping testing of sex specific effects."
                touch ${prefix}NOT_DONE.sex_spec.gz
                continue
            fi

            base=${prefix}_\${p}.regenie.gz
            zcat \$base | extract_sex_variants.awk -v threshold=${sex_logp} > \${p}.sex_variants
            nvars=\$(wc -l < \${p}.sex_variants)
            echo "\$nvars variants will be tested for sex specific effects"

            if [[ \$nvars -lt 1 ]]; then
                write_empty_sex_spec_header.sh ${is_binary} ${prefix}.\${p}.sex_spec.gz
                continue
            fi

            for s in males females; do
                echo "running \$s analysis"
                echo "\$(wc -l \$s) individuals in file"
                head \$s

                regenie \\
                    --step 2 \\
                    ${test_cmd} \\
                    ${binary_flag} \\
                    ${geno_flags} \\
                    --keep \$s \\
                    --extract \${p}.sex_variants \\
                    --covarFile ${cov_pheno} \\
                    --covarColList \$sex_covars \\
                    --phenoFile ${cov_pheno} \\
                    --phenoColList \$p \\
                    --pred ${pred} \\
                    --bsize ${bsize} \\
                    --threads \$n_cpu \\
                    --gz \\
                    --out ${prefix}.sex_spec.\$s \\
                    ${options}
            done

            # regenie writes *.sex_spec.<males|females>_<pheno>.regenie.gz —
            # merge_sex_specific.py expects *.sex_spec.<males|females>_<pheno>.gz
            for f in ${prefix}.sex_spec.*_\${p}.regenie.gz; do
                mv \$f \${f%.regenie.gz}.gz
            done

            merge_sex_specific.py \$p ${prefix} ${prefix}_\${p}.regenie.gz
        done
    else
        # Sex-specific disabled: emit a single placeholder so the output channel
        # stays non-empty and GATHER can detect no real shards exist.
        touch ${prefix}NOT_DONE.sex_spec.gz
    fi
    """

    stub:
    geno_tag = geno.getName()
    prefix   = pred.getName().replaceAll(/_pred\.list$/, '') + '.' + geno_tag
    """
    touch ${prefix}.log
    ${phenolist.collect { "touch ${prefix}_${it}.regenie.gz" }.join('\n    ')}
    touch ${prefix}NOT_DONE.sex_spec.gz
    """
}
