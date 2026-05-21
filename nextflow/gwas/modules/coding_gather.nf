// Global coding-gather: concatenates every pheno's <pheno>_coding.txt into a
// single bgzipped coding_variants.txt.gz. Mirrors wdl/gwas/regenie.wdl
// (coding_gather task).

process CODING_GATHER {
    tag "coding_gather"
    container params.coding_gather_docker

    cpus   1
    memory 2.GB
    disk   20.GB, type: 'pd-standard'

    publishDir "${params.outdir}", mode: 'copy', pattern: 'coding_variants.txt.gz'

    input:
    path(coding_files, stageAs: 'coding/*')

    output:
    path("coding_variants.txt.gz"), emit: main

    script:
    """
    set -euxo pipefail
    files=( coding/* )
    cat <(head -n1 \${files[0]}) <(awk 'FNR>1' "\${files[@]}") | bgzip > coding_variants.txt.gz
    """

    stub:
    """
    touch coding_variants.txt.gz
    """
}
