#!/bin/bash
# Combine GWAS results across ALL phenotypes (bgen autosomes + pgen X/Y/MT),
# filter to mlogp >= -log10(pval_max), and prepend a phenotype column.
# Streams directly from GCS, no per-phenotype files written.
#
# Usage: ./combine_all_filtered.sh <output_file.gz> [parallel_jobs] [pval_threshold]
#   output_file.gz  local path or gs:// path (single combined gzipped TSV)
#   parallel_jobs   phenotypes streamed concurrently (default 16)
#   pval_threshold  keep rows with pval < this (default 1e-4)

set -euo pipefail

BGEN_DIR="gs://r14-data/chip_analysis/additive_chip_gwas/nf-outputs-bgen/munged/munged"
PGEN_DIR="gs://r14-data/chip_analysis/additive_chip_gwas/nf-outputs-pgen/munged/munged"

OUT_FILE="${1:?usage: $0 <output_file.gz> [parallel_jobs] [pval_threshold]}"
JOBS="${2:-16}"
PVAL_MAX="${3:-1e-4}"

# Use mlogp (column 6) for filtering: more robust against tiny/NA pvals than parsing pval text.
MLOGP_MIN=$(awk -v p="$PVAL_MAX" 'BEGIN{ printf "%.10g", -log(p)/log(10) }')
echo "Filter: pval < ${PVAL_MAX}  =>  mlogp > ${MLOGP_MIN}" >&2

TMP_DIR="$(mktemp -d -t combine_filtered.XXXXXX)"
trap 'rm -rf "$TMP_DIR"' EXIT

mapfile -t PHENOS < <(
    comm -12 \
        <(gsutil ls "${BGEN_DIR}/" | grep -oP '[^/]+\.gz$' | sed 's/\.gz$//' | sort -u) \
        <(gsutil ls "${PGEN_DIR}/" | grep -oP '[^/]+\.gz$' | sed 's/\.gz$//' | sort -u)
)

echo "Filtering ${#PHENOS[@]} phenotypes (parallel=${JOBS}) -> ${OUT_FILE}" >&2

process_one() {
    local pheno="$1"
    local out_chunk="${TMP_DIR}/${pheno}.tsv"
    # Concatenated gzip streams decompress as one stream; awk drops headers, filters,
    # and remaps X->23, Y->24 (MT already 25). No per-chunk sort -- one global sort below.
    { gsutil cat "${BGEN_DIR}/${pheno}.gz"
      gsutil cat "${PGEN_DIR}/${pheno}.gz"
    } | zcat \
      | awk -v pheno="$pheno" -v mmin="$MLOGP_MIN" 'BEGIN{FS=OFS="\t"}
            /^#/ { next }
            $6 != "NA" && $6+0 > mmin+0 {
                if ($1=="X") $1=23; else if ($1=="Y") $1=24;
                print pheno, $0
            }' \
      > "$out_chunk"
}

export -f process_one
export BGEN_DIR PGEN_DIR TMP_DIR MLOGP_MIN

printf '%s\n' "${PHENOS[@]}" \
    | xargs -n1 -P"$JOBS" -I{} bash -c 'process_one "$@"' _ {}

echo "Merging chunks -> ${OUT_FILE}" >&2

merged_stream() {
    printf 'pheno\t#chrom\tpos\tref\talt\tpval\tmlogp\tbeta\tsebeta\taf_alt\taf_alt_cases\taf_alt_controls\n'
    # Use find+cat so the arg list scales beyond shell glob limits.
    # Global sort by (chrom, pos, ref, alt) with pheno as deterministic tiebreaker.
    find "$TMP_DIR" -maxdepth 1 -name '*.tsv' -print0 \
        | xargs -0 cat \
        | LC_ALL=C sort -T "$TMP_DIR" -t$'\t' -k2,2n -k3,3n -k4,4 -k5,5 -k1,1
}

if [[ "$OUT_FILE" == gs://* ]]; then
    merged_stream | gzip | gsutil -q cp - "$OUT_FILE"
else
    merged_stream | gzip > "${OUT_FILE}.tmp" && mv "${OUT_FILE}.tmp" "$OUT_FILE"
fi

echo "Wrote ${OUT_FILE}" >&2
