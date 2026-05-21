// Per-pheno FinnGen annotation wrapper. Calls annotate_summary.py to produce
// pheweb-compatible <pheno>_summary.txt and <pheno>_coding.txt. Mirrors
// wdl/gwas/regenie_sub.wdl (summary task).

process SUMMARY {
    tag "${pheno}"
    container params.summary_docker

    cpus   1
    memory 2.GB
    disk   200.GB, type: 'pd-standard'

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: '*_summary.txt'
    publishDir "${params.outdir}/coding",  mode: 'copy', pattern: '*_coding.txt'

    input:
    tuple val(pheno),
          path(pheweb),
          path(finngen_annotation),
          path(finngen_annotation_tbi)

    output:
    tuple val(pheno), path("*_summary.txt"), path("*_coding.txt"), emit: main

    script:
    """
    set -euxo pipefail
    annotate_summary.py ${pheweb} ${finngen_annotation} ${params.summary_pval_thresh} ${params.coding_pval_thresh}
    """

    stub:
    """
    base=\$(basename ${pheweb} .gz)
    touch \${base}_summary.txt \${base}_coding.txt
    """
}
