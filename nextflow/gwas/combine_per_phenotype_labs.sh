#!/bin/bash
# Combine bgen (autosomes) + pgen (X/Y/MT) munged GWAS results per phenotype.
# Output is one gzipped file per phenotype written to <output_dir>.
#
# Usage: ./combine_per_phenotype.sh <output_dir> [parallel_jobs]
#   output_dir     local path or gs:// prefix
#   parallel_jobs  number of phenotypes to process concurrently (default 8)

set -euo pipefail

BGEN_DIR="gs://r14-data/chip_analysis/labs_quant_gwas/nf-outputs-bgen/munged/munged"
PGEN_DIR="gs://r14-data/chip_analysis/labs_quant_gwas/nf-outputs-pgen/munged/munged"

OUT_DIR="${1:?usage: $0 <output_dir> [parallel_jobs]}"
JOBS="${2:-8}"

OUT_DIR="${OUT_DIR%/}"
if [[ "$OUT_DIR" != gs://* ]]; then
    mkdir -p "$OUT_DIR"
fi

# Intersect phenotype lists (data files only, not .tbi indexes).
mapfile -t PHENOS < <(
    comm -12 \
        <(gsutil ls "${BGEN_DIR}/" | grep -oP '[^/]+\.gz$' | sed 's/\.gz$//' | sort -u) \
        <(gsutil ls "${PGEN_DIR}/" | grep -oP '[^/]+\.gz$' | sed 's/\.gz$//' | sort -u)
)

echo "Combining ${#PHENOS[@]} phenotypes (parallel=${JOBS}) into ${OUT_DIR}/" >&2

HEADER=$'#chrom\tpos\tref\talt\tpval\tmlogp\tbeta\tsebeta\taf_alt\taf_alt_cases\taf_alt_controls'

combine_one() {
    local pheno="$1"
    local src_bgen="${BGEN_DIR}/${pheno}.gz"
    local src_pgen="${PGEN_DIR}/${pheno}.gz"
    local dst="${OUT_DIR}/${pheno}.gz"

    # Remap pgen X->23, Y->24 (MT already 25), then global sort by (chrom, pos, ref, alt).
    # Multi-allelic sites in the bgen source aren't allele-ordered, so we sort the whole
    # body — not just the pgen block. Header is fixed-schema, prepended after sort.
    build_stream() {
        printf '%s\n' "$HEADER"
        { gsutil cat "$src_bgen" | zcat | tail -n +2
          gsutil cat "$src_pgen" | zcat | tail -n +2 \
              | awk 'BEGIN{FS=OFS="\t"} {if($1=="X")$1=23; else if($1=="Y")$1=24; print}'
        } | LC_ALL=C sort -t$'\t' -k1,1n -k2,2n -k3,3 -k4,4
    }

    if [[ "$dst" == gs://* ]]; then
        build_stream | gzip | gsutil -q cp - "$dst"
    else
        build_stream | gzip > "${dst}.tmp" && mv "${dst}.tmp" "$dst"
    fi
    echo "done: $pheno" >&2
}

export -f combine_one
export BGEN_DIR PGEN_DIR OUT_DIR HEADER

printf '%s\n' "${PHENOS[@]}" \
    | xargs -n1 -P"$JOBS" -I{} bash -c 'combine_one "$@"' _ {}

echo "All ${#PHENOS[@]} phenotypes written to ${OUT_DIR}/" >&2
