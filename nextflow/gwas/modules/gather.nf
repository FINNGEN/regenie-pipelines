// Per-pheno gather: concatenates regenie shards across genotype chunks,
// sorts them, bgzips, generates a pheweb-format munged file, and runs qqplot.
// Also concatenates sex-specific shards if any real ones exist. Mirrors
// wdl/gwas/regenie_sub.wdl (gather task).

process GATHER {
    tag "${pheno}"
    container params.gather_docker

    cpus   1
    memory 8.GB
    disk   200.GB, type: 'pd-standard'

    publishDir "${params.outdir}/regenie",    mode: 'copy', pattern: 'regenie/*.gz*'
    publishDir "${params.outdir}/munged",     mode: 'copy', pattern: 'munged/*.gz*'
    // Disambiguate per-pheno on publish: every task writes the same in-workdir
    // names (qq_out / qq_err and the skip_qq NOTDONE placeholders), so many
    // tasks racing to copy to the same gs:// destination trips publishDir.
    publishDir "${params.outdir}/qq",         mode: 'copy', pattern: '*.png',           saveAs: { f -> f.startsWith('NOTDONE') ? "${pheno}.${f}" : f }
    publishDir "${params.outdir}/qq",         mode: 'copy', pattern: '*qquantiles.txt', saveAs: { f -> f.startsWith('NOTDONE') ? "${pheno}.${f}" : f }
    publishDir "${params.outdir}/qq",         mode: 'copy', pattern: 'qq_out',          saveAs: { "${pheno}.qq_out" }
    publishDir "${params.outdir}/qq",         mode: 'copy', pattern: 'qq_err',          saveAs: { "${pheno}.qq_err" }

    input:
    tuple val(pheno),
          path(regenie_shards, stageAs: 'shards/*'),
          path(sex_shards,     stageAs: 'sex/*')

    output:
    tuple val(pheno),
          path("regenie/${pheno}.regenie.gz"),
          path("regenie/${pheno}.regenie.sex_diff.gz"),
          path("regenie/${pheno}.regenie.sex_diff.gz.tbi"),
          path("munged/${pheno}.gz"),
          path("munged/${pheno}.gz.tbi"),
          path("qq_out"),
          path("qq_err"),
          path("*.png"),
          path("*qquantiles.txt"), emit: main

    script:
    def is_binary = params.is_binary
    def sort_flag = params.gather_sort_flag
    def skip_qq   = params.skip_qqplot as boolean
    def munger    = is_binary ? 'munge_regenie_binary.awk' : 'munge_regenie_quant.awk'
    """
    set -euxo pipefail
    mkdir -p regenie munged

    # 1) Concatenate regenie shards (sorted)
    shards=( shards/* )
    echo -e "\$(date)\tconcatenating \${#shards[@]} regenie shards into regenie/${pheno}.regenie.gz"
    cat \\
        <(zcat \${shards[0]} | head -1) \\
        <(for f in "\${shards[@]}"; do zcat \$f | tail -n+2; done | sort ${sort_flag}) \\
        | bgzip > regenie/${pheno}.regenie.gz

    # 2) Concatenate sex-specific shards — only if any real (non-NOT_DONE) shards exist.
    real_sex=()
    if [[ -d sex ]]; then
        for f in sex/*; do
            [[ -e \$f ]] || continue
            base=\$(basename \$f)
            if [[ "\$base" == *NOT_DONE* ]]; then
                continue
            fi
            real_sex+=("\$f")
        done
    fi

    if [[ \${#real_sex[@]} -ge 1 ]]; then
        echo -e "\$(date)\tconcatenating \${#real_sex[@]} sex-specific shards"
        cat \\
            <(zcat \${real_sex[0]} | awk '{ print "#"\$0; exit 0 }') \\
            <(for f in "\${real_sex[@]}"; do zcat \$f | tail -n+2; done | sort ${sort_flag}) \\
            | tr ' ' '\\t' | bgzip > regenie/${pheno}.regenie.sex_diff.gz
        tabix -s 1 -b 2 -e 2 regenie/${pheno}.regenie.sex_diff.gz
    else
        touch regenie/${pheno}.regenie.sex_diff.gz
        touch regenie/${pheno}.regenie.sex_diff.gz.tbi
    fi

    # 3) Munge to pheweb format
    echo -e "\$(date)\tmunging regenie/${pheno}.regenie.gz → munged/${pheno}.gz"
    zcat regenie/${pheno}.regenie.gz | ${munger} | bgzip > munged/${pheno}.gz

    # 4) QQ-plot (skippable for chip pgen variant)
    if [[ "${skip_qq}" == "true" ]]; then
        : > qq_out
        : > qq_err
        touch NOTDONE.png NOTDONEqquantiles.txt
    else
        echo -e "\$(date)\tpreparing qqplot input"
        zcat munged/${pheno}.gz | cut -f1,2,5 \\
            | awk 'BEGIN{ FS=OFS="\\t" } NR==1{ print \$0,"SNP" } NR>1{ print \$0,NR }' > ${pheno}
        head ${pheno}
        echo -e "\$(date)\trunning qqplot.R"
        qqplot.R --file ${pheno} --chrcol "#chrom" --bp_col "pos" --pval_col "pval" --loglog_pval 10 > qq_out 2> qq_err
    fi

    # 5) Tabix the munged file
    echo -e "\$(date)\ttabixing munged/${pheno}.gz"
    tabix -s 1 -b 2 -e 2 munged/${pheno}.gz
    echo -e "\$(date)\tdone"
    """

    stub:
    """
    mkdir -p regenie munged
    touch regenie/${pheno}.regenie.gz
    touch regenie/${pheno}.regenie.sex_diff.gz
    touch regenie/${pheno}.regenie.sex_diff.gz.tbi
    touch munged/${pheno}.gz
    touch munged/${pheno}.gz.tbi
    touch qq_out qq_err
    touch STUB.png STUBqquantiles.txt
    """
}
