#!/usr/bin/env bash
# Bash port of filter_hits_regions.py -- same operations, for local testing without python, shaped
# to paste directly into extract_cond_regions's WDL command block (hardcoded inputs instead of a
# CLI, matching how regenie_conditional_bash was inlined). See filter_hits_regions.py for the
# authoritative implementation; keep the two in sync.
#
# For each region (a "chrom start end" bed-like file, built below from the pheno's finemap regions
# plus HLA/PAR), finds the single most significant variant (by mlogp) within that region from
# --sumstats, restricted to chromosomes in CHROMS. Writes one combined hits file (pheno, chrom,
# region, variant, mlogp) plus one file per chromosome.
set -euo pipefail

# Hardcoded, WDL-task-shaped inputs -- no CLI parsing, matching how extract_cond_regions passes
# these in (root paths + PHENO substitution). Edit these directly to test a different pheno.
PHENO="T2D"
SUMSTATS="$HOME/r14/regenie/release/summary_stats/${PHENO}.gz"
REGION_ROOT="$HOME/r14/finemap/release/beds/${PHENO}.bed"
PVAL_THRESHOLD=7.3
CHR_COL="#chrom"
POS_COL="pos"
REF_COL="ref"
ALT_COL="alt"
MLOGP_COL="mlogp"
CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23)
OUT="."

[[ -f "$SUMSTATS" ]] || { echo "sumstats file not found: $SUMSTATS" >&2; exit 1; }
[[ -f "$REGION_ROOT" ]] || { echo "region bed not found: $REGION_ROOT" >&2; exit 1; }

mkdir -p "$OUT"

# same as extract_cond_regions's own region-building: start from the pheno's finemap regions,
# then append extra fixed regions -- HLA, and PAR (pseudoautosomal, chrX).
REGIONS="tmp.bed"
cat "$REGION_ROOT" > "$REGIONS"
echo -e "6\t29000000\t34000000" >> "$REGIONS"
# TODO: PAR region coordinates -- need exact chrom/start/end before adding.
echo "$(IFS=' '; echo "${CHROMS[*]}") chromosomes included."
echo "> ${PHENO} <"

REGION_TMP="${OUT}/${PHENO}_subset_regions.txt.region.bed"
CHROM_LIST=$(IFS=,; echo "${CHROMS[*]}")

# keep only regions whose chrom (column 1, already a bare digit string in the finemap bed) is in
# CHROMS -- written out UNMODIFIED, same as filter_hits_regions.py's o.write(line)
awk -v chrom_list="${CHROM_LIST}" '
  BEGIN { n = split(chrom_list, arr, ","); for (i=1;i<=n;i++) valid[arr[i]] = 1 }
  ($1 in valid)
' "$REGIONS" > "$REGION_TMP"

TOT_REGIONS=$(wc -l < "$REGION_TMP")
[[ "$TOT_REGIONS" -gt 0 ]] && echo "${TOT_REGIONS} regions to be filtered." || echo "WARNING: ${REGIONS} empty. No hits will be returned!" >&2

echo -n "Processing data..."

# column indices, read once from the sumstats header (tabix region syntax is 1-based/inclusive on
# both ends -- same convention the sumstats pos column already uses, so no coordinate-system
# mismatch to worry about, unlike e.g. BED-based tools)
HEADER=$(tabix -H "$SUMSTATS")
read -r P M R A <<< "$(awk -F'\t' -v pos_col="$POS_COL" -v mlogp_col="$MLOGP_COL" -v ref_col="$REF_COL" -v alt_col="$ALT_COL" '
  { for(i=1;i<=NF;i++) idx[$i]=i; print idx[pos_col], idx[mlogp_col], idx[ref_col], idx[alt_col] }
' <<< "$HEADER")"

OUT_FILE="${OUT}/${PHENO}_sig_hits.txt"
rm -f "${OUT}/${PHENO}_sig_hits_"*.txt "$OUT_FILE"
# always create OUT_FILE, even with zero hits -- matches python's unconditional `with open(...)`.
# the WDL declares this as a required (non-optional) File output, so a shard with no hits above
# threshold must still produce an (empty) file or the task fails.
touch "$OUT_FILE"

while read -r chrom start end; do
  # sort's output can be large enough that `head -1` closes the pipe before sort finishes writing,
  # killing sort with SIGPIPE (rc 141) -- pipefail then propagates that and set -e kills the whole
  # script. Drain the rest after head takes its line so sort always exits cleanly.
  top=$(tabix "$SUMSTATS" "${chrom}:${start}-${end}" \
    | awk -F'\t' -v m="$M" -v t="$PVAL_THRESHOLD" '$m>t' \
    | sort -t$'\t' -k"${M},${M}rn" | { head -1; cat >/dev/null; })
  [[ -z "$top" ]] && continue

  IFS=$'\t' read -r -a f <<< "$top"
  mlogp="${f[$((M-1))]}" pos="${f[$((P-1))]}" ref="${f[$((R-1))]}" alt="${f[$((A-1))]}"
  # matches the python script inv_sub_dict: chrom 23 is written back out as X in variant IDs
  # (sumstats/regions use bare "23", but variant-id convention uses "chrX")
  label="$chrom"; [[ "$chrom" == "23" ]] && label="X"
  variant="chr${label}_${pos}_${ref}_${alt}"
  region="${chrom}:${start}-${end}"
  line="${PHENO}\t${chrom}\t${region}\t${variant}\t${mlogp}"
  echo -e "$line" >> "$OUT_FILE"
  echo -e "$line" >> "${OUT}/${PHENO}_sig_hits_${chrom}.txt"
done < "$REGION_TMP"

echo "done."
echo "dumping results to ${OUT_FILE}"
