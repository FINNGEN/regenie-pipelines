#!/bin/bash
# Combine bgen (autosomes) + pgen (X/Y/MT) GWAS results for the three
# quant phenotypes whose pgen run mislabels chr25 (MT) and chrY rows as
# CHROM=23 in the regenie output (a regenie quirk: without --bt, non-X
# non-autosomal codes collapse to 23). The variant ID column still carries
# the true chromosome, so we re-munge from the gathered raw regenie file
# (.../regenie/regenie/<pheno>.regenie.gz) using ID-derived CHROM.
#
# Output is 11-column pheweb format with NA padding for the case/control
# AF columns that quants don't have, so the file lines up with the binary
# combined files for downstream concatenation:
#   #chrom pos ref alt pval mlogp beta sebeta af_alt af_alt_cases af_alt_controls
#
# Usage: ./combine_per_phenotype_subset.sh <output_dir> [parallel_jobs]

set -euo pipefail

BGEN_DIR="gs://r14-data/chip_analysis/additive_chip_gwas/nf-outputs-bgen/munged/munged"
PGEN_REGENIE_DIR="gs://r14-data/chip_analysis/additive_chip_gwas/nf-outputs-pgen/regenie/regenie"

OUT_DIR="${1:?usage: $0 <output_dir> [parallel_jobs]}"
JOBS="${2:-3}"

OUT_DIR="${OUT_DIR%/}"
if [[ "$OUT_DIR" != gs://* ]]; then
    mkdir -p "$OUT_DIR"
fi

PHENOS=(BMI_IRN HEIGHT_IRN WEIGHT_IRN)

echo "Combining ${#PHENOS[@]} phenotypes (parallel=${JOBS}) into ${OUT_DIR}/" >&2

HEADER=$'#chrom\tpos\tref\talt\tpval\tmlogp\tbeta\tsebeta\taf_alt\taf_alt_cases\taf_alt_controls'

combine_one() {
    local pheno="$1"
    local src_bgen="${BGEN_DIR}/${pheno}.gz"
    local src_pgen_regenie="${PGEN_REGENIE_DIR}/${pheno}.regenie.gz"
    local dst="${OUT_DIR}/${pheno}.gz"

    # Re-munge pgen-side from raw regenie output, taking CHROM from the
    # variant ID (chrX→23, chrY→24, chrM/chrMT/chr25→25, chrN→N). Drop
    # rows with LOG10P=NA, mirroring the standard quant munger behaviour
    # for failed tests (defensive — regenie quant rarely emits NA but
    # firth-se can).
    build_stream() {
        printf '%s\n' "$HEADER"
        { gsutil cat "$src_bgen" | zcat | tail -n +2 \
              | awk 'BEGIN{FS=OFS="\t"} {print $0, "NA", "NA"}'
          gsutil cat "$src_pgen_regenie" | zcat \
              | awk 'BEGIN{OFS="\t"}
                     NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
                     $h["LOG10P"] == "NA" { next }
                     {
                         split($h["ID"], a, "_"); c=a[1]; sub(/^chr/, "", c)
                         if (c=="X") c=23
                         else if (c=="Y") c=24
                         else if (c=="M" || c=="MT") c=25
                         print c, $h["GENPOS"], $h["ALLELE0"], $h["ALLELE1"],
                               10^-$h["LOG10P"], $h["LOG10P"],
                               $h["BETA"], $h["SE"], $h["A1FREQ"], "NA", "NA"
                     }'
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
export BGEN_DIR PGEN_REGENIE_DIR OUT_DIR HEADER

printf '%s\n' "${PHENOS[@]}" \
    | xargs -n1 -P"$JOBS" -I{} bash -c 'combine_one "$@"' _ {}

echo "All ${#PHENOS[@]} phenotypes written to ${OUT_DIR}/" >&2
