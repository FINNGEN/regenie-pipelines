#!/bin/bash
# Bash port of regenie_conditional.py -- same CLI, same operations, for local testing without pandas/python.
# See regenie_conditional.py for the authoritative implementation; keep the two in sync.
set -euo pipefail

REGENIE_COVARIATES="SEX_IMPUTED,AGE_AT_DEATH_OR_END_OF_FOLLOWUP,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,IS_FINNGEN2_CHIP,BATCH_DS1_BOTNIA_Dgi_norm,BATCH_DS10_FINRISK_Palotie_norm,BATCH_DS11_FINRISK_PredictCVD_COROGENE_Tarto_norm,BATCH_DS12_FINRISK_Summit_norm,BATCH_DS13_FINRISK_Bf_norm,BATCH_DS14_GENERISK_norm,BATCH_DS15_H2000_Broad_norm,BATCH_DS16_H2000_Fimm_norm,BATCH_DS17_H2000_Genmets_norm_relift,BATCH_DS18_MIGRAINE_1_norm_relift,BATCH_DS19_MIGRAINE_2_norm,BATCH_DS2_BOTNIA_T2dgo_norm,BATCH_DS20_SUPER_1_norm_relift,BATCH_DS21_SUPER_2_norm_relift,BATCH_DS22_TWINS_1_norm,BATCH_DS23_TWINS_2_norm_nosymmetric,BATCH_DS24_SUPER_3_norm,BATCH_DS25_BOTNIA_Regeneron_norm,BATCH_DS26_DIREVA_norm,BATCH_DS27_NFBC66_norm,BATCH_DS28_NFBC86_norm,BATCH_DS3_COROGENE_Sanger_norm,BATCH_DS4_FINRISK_Corogene_norm,BATCH_DS5_FINRISK_Engage_norm,BATCH_DS6_FINRISK_FR02_Broad_norm_relift,BATCH_DS7_FINRISK_FR12_norm,BATCH_DS8_FINRISK_Finpcga_norm,BATCH_DS9_FINRISK_Mrpred_norm"

PVAL_THRESHOLD=7
COVARIATES="$REGENIE_COVARIATES"
REGENIE_PARAMS=" --bt --firth --approx --bsize 200 --ref-first"
REGENIE_CMD="regenie"
FORCE=0
MAX_STEPS=10
CHR_COL="#chrom"
POS_COL="pos"
REF_COL="ref"
ALT_COL="alt"
MLOGP_COL="mlogp"
BETA_COL="beta"
SEBETA_COL="sebeta"
THREADS=$(nproc)
SAMPLE_FILE=""
LOCUS_ARG=""
REGION_ARG=""
LOCUS_LIST=""
PHENO=""
OUT=""
PHENO_FILE="/home/pete/r14/pheno/R14_COV_PHENO_V0.FID.txt.gz"
BGEN=""
SUMSTATS=""
NULL_FILE=""

usage() {
  cat >&2 <<'EOF'
Usage: regenie_conditional.sh --pheno P --out OUT --bgen B --sumstats S --null-file N
                               (--locus-region LOCUS REGION | --locus-list FILE) [options]
  --pval-threshold FLOAT     threshold limit (-log10(p)), or a raw p-value (<1) (default 7)
  --pheno-file FILE          default: /home/pete/r14/pheno/R14_COV_PHENO_V0.FID.txt.gz
  --covariates LIST          comma-separated covariate list (default: full FinnGen list)
  --sample-file FILE         bgen sample file (auto-detected next to --bgen if omitted)
  --regenie-params STR       extra regenie params (default: " --bt --firth --approx --bsize 200 --ref-first")
  --force                    force re-run of already-completed steps
  --max-steps INT            default 10
  --chr-col/--pos-col/--ref-col/--alt-col/--mlogp-col/--beta-col/--sebeta-col STR
  --threads INT               default: nproc
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pval-threshold) PVAL_THRESHOLD="$2"; shift 2 ;;
    --pheno) PHENO="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --covariates) COVARIATES="$2"; shift 2 ;;
    --pheno-file) PHENO_FILE="$2"; shift 2 ;;
    --bgen) BGEN="$2"; shift 2 ;;
    --sample-file) SAMPLE_FILE="$2"; shift 2 ;;
    --sumstats) SUMSTATS="$2"; shift 2 ;;
    --regenie-params) REGENIE_PARAMS="$2"; shift 2 ;;
    --null-file) NULL_FILE="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    --max-steps) MAX_STEPS="$2"; shift 2 ;;
    --chr_col|--chr-col) CHR_COL="$2"; shift 2 ;;
    --pos_col|--pos-col) POS_COL="$2"; shift 2 ;;
    --ref_col|--ref-col) REF_COL="$2"; shift 2 ;;
    --alt_col|--alt-col) ALT_COL="$2"; shift 2 ;;
    --mlogp_col|--mlogp-col) MLOGP_COL="$2"; shift 2 ;;
    --beta_col|--beta-col) BETA_COL="$2"; shift 2 ;;
    --sebeta_col|--sebeta-col) SEBETA_COL="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --locus-region) LOCUS_ARG="$2"; REGION_ARG="$3"; shift 3 ;;
    --locus-list) LOCUS_LIST="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown argument: $1" >&2; usage ;;
  esac
done

[[ -n "$PHENO" ]] || { echo "--pheno is required" >&2; usage; }
[[ -n "$OUT" ]] || { echo "--out is required" >&2; usage; }
[[ -n "$PHENO_FILE" && -f "$PHENO_FILE" ]] || { echo "--pheno-file must exist (default: /home/pete/r14/pheno/R14_COV_PHENO_V0.FID.txt.gz)" >&2; usage; }
[[ -n "$BGEN" && -f "$BGEN" ]] || { echo "--bgen is required and must exist" >&2; usage; }
[[ -n "$SUMSTATS" && -f "$SUMSTATS" ]] || { echo "--sumstats is required and must exist" >&2; usage; }
[[ -n "$NULL_FILE" && -f "$NULL_FILE" ]] || { echo "--null-file is required and must exist" >&2; usage; }
if [[ -n "$LOCUS_ARG" && -n "$LOCUS_LIST" ]] || [[ -z "$LOCUS_ARG" && -z "$LOCUS_LIST" ]]; then
  echo "exactly one of --locus-region or --locus-list is required" >&2; usage
fi
[[ -z "$LOCUS_LIST" || -f "$LOCUS_LIST" ]] || { echo "--locus-list file does not exist" >&2; usage; }
[[ -z "$SAMPLE_FILE" || -f "$SAMPLE_FILE" ]] || { echo "--sample-file does not exist" >&2; usage; }

pretty_print() {
  local string="$1" l="${2:-30}"
  local half=$(( l - ${#string}/2 )); (( half < 0 )) && half=0
  local dashes; dashes=$(printf '%*s' "$half" '' | tr ' ' '-')
  echo "${dashes}> ${string} <${dashes}"
}

map_vals_to_string() {
  awk -v b="$1" -v s="$2" -v m="$3" 'BEGIN{printf "{'"'"'beta'"'"': %.2f, '"'"'sebeta'"'"': %.2f, '"'"'mlogp'"'"': %.2f}", b, s, m}'
}

# region[0] can be "locus" or "region" depending on argument order -- whichever looks like CHR:START-END is the region.
check_region() {
  local locus="$1" region="$2"
  if [[ "$locus" == *:* && "$locus" == *-* ]]; then
    local tmp="$locus"; locus="$region"; region="$tmp"
  fi
  CHECKED_LOCUS="$locus"
  CHECKED_REGION_FLAG=" --range $region "
}

OUT_DIR=$(dirname "$OUT")
BASENAME=$(basename "$OUT")
TMP_DIR="$OUT_DIR/tmp"
mkdir -p "$OUT_DIR" "$TMP_DIR"

# automatically look for sample file if not explicitly passed
if [[ -z "$SAMPLE_FILE" ]]; then
  for candidate in "${BGEN%.bgen}.bgen.sample" "${BGEN%.bgen}.sample"; do
    [[ -f "$candidate" ]] && { SAMPLE_FILE="$candidate"; echo "WARNING: Using $SAMPLE_FILE as sample file." >&2; break; }
  done
fi

# figure out whether threshold is a p-value or -log10(p)
PVAL_THRESHOLD=$(awk -v v="$PVAL_THRESHOLD" 'BEGIN{printf "%.10g", (v<1) ? -log(v)/log(10) : v}')
pretty_print "MLOGP THRESHOLD: $PVAL_THRESHOLD"

# ---- sumstats lookup cache: snpid -> beta, sebeta, mlogp (from the ORIGINAL, unconditioned sumstats) ----
PVAL_CACHE="$TMP_DIR/${PHENO}_pvals.tsv"
if [[ ! -f "$PVAL_CACHE" || "$FORCE" -eq 1 ]]; then
  echo "reading original pvals.." >&2
  zcat -f "$SUMSTATS" | awk -F'\t' -v chr_col="$CHR_COL" -v pos_col="$POS_COL" -v ref_col="$REF_COL" -v alt_col="$ALT_COL" -v mlogp_col="$MLOGP_COL" -v beta_col="$BETA_COL" -v sebeta_col="$SEBETA_COL" '
    NR==1 {
      for (i=1;i<=NF;i++) idx[$i]=i
      c=idx[chr_col]; p=idx[pos_col]; r=idx[ref_col]; a=idx[alt_col]; m=idx[mlogp_col]; b=idx[beta_col]; s=idx[sebeta_col]
      next
    }
    { print "chr" $c "_" $p "_" $r "_" $a "\t" $b "\t" $s "\t" $m }
  ' > "$PVAL_CACHE"
  echo "done." >&2
fi

# get_sum_dict_data equivalent: prints "beta\tsebeta\tmlogp" for a variant, or "0\t0\t0" if not found
# (matches the python script's defaultdict(lambda: defaultdict(lambda: "0")) fallback)
get_sum_dict_data() {
  local variant="$1"
  local line
  line=$(awk -F'\t' -v v="$variant" '$1==v{print $2"\t"$3"\t"$4; exit}' "$PVAL_CACHE")
  [[ -n "$line" ]] && echo "$line" || echo -e "0\t0\t0"
}

# ---- filter_pheno: extract FID/IID + pheno + covariate columns into a small tmp pheno file ----
FILTERED_PHENO_FILE="$TMP_DIR/${PHENO}_pheno.tmp"
if [[ ! -f "$FILTERED_PHENO_FILE" || "$FORCE" -eq 1 ]]; then
  echo "Creating new pheno file..." >&2
  header=$(zcat -f "$PHENO_FILE" | { head -1; cat >/dev/null; })
  cols=$(awk -F'\t' -v want="${PHENO},${COVARIATES}" -v h="$header" '
    BEGIN {
      n = split(want, w, ",")
      split(h, f, "\t")
      for (i in f) idx[f[i]] = i
      out = "1,2"
      for (j=1; j<=n; j++) {
        if (!(w[j] in idx)) { print "column " w[j] " not found in pheno file" > "/dev/stderr"; exit 1 }
        out = out "," idx[w[j]]
      }
      print out
    }')
  zcat -f "$PHENO_FILE" | cut -f"$cols" > "$FILTERED_PHENO_FILE"
  echo "done." >&2
fi

# ---- single regenie conditional run; sets REGENIE_OUT_FILE and REGENIE_RET ----
regenie_run() {
  local step="$1" locus="$2" region_flag="$3" log_file="$4"

  local regenie_file="$TMP_DIR/${BASENAME}_${PHENO}.regenie"
  local out_file="$OUT_DIR/${BASENAME}_${PHENO}_${locus}_${step}.conditional"
  pretty_print "VARIANT:${CONDITION_LIST[-1]}"
  echo "generating $out_file ..." >&2

  if [[ ! -f "$out_file" || "$FORCE" -eq 1 ]]; then
    FORCE=1

    local pred_file="$TMP_DIR/${BASENAME}_${PHENO}.pred"
    printf "%s\t%s\n" "$PHENO" "$NULL_FILE" > "$pred_file"

    local tmp_variant="$TMP_DIR/${BASENAME}_variant.tmp"
    printf "%s\n" "${CONDITION_LIST[@]}" > "$tmp_variant"

    local sample_cmd=""
    [[ -n "$SAMPLE_FILE" && -f "$SAMPLE_FILE" ]] && sample_cmd=" --sample $SAMPLE_FILE"

    local cmd="$REGENIE_CMD --step 2 $REGENIE_PARAMS --bgen $BGEN${sample_cmd} --out $TMP_DIR/$BASENAME --pred $pred_file --phenoFile $FILTERED_PHENO_FILE --phenoCol $PHENO --condition-list $tmp_variant $region_flag --covarFile $FILTERED_PHENO_FILE --covarColList $COVARIATES --threads $THREADS"
    echo "$cmd" >&2

    local start=$SECONDS ret=0
    bash -c "$cmd" >> "$log_file" || ret=$?
    echo "Script ran in $((SECONDS-start)) seconds with $THREADS cpus." >&2

    [[ "$ret" -eq 0 ]] && mv "$regenie_file" "$out_file"
    REGENIE_OUT_FILE="$out_file"
    REGENIE_RET="$ret"
  else
    echo "file already exists" >&2
    REGENIE_OUT_FILE="$out_file"
    REGENIE_RET=0
  fi
}

# ---- check_hit: find the top (max LOG10P) row in a regenie output file; sets NEXT_STEP and HIT_* ----
check_hit() {
  local out_file="$1" step="$2" threshold="$3"
  local result
  result=$(awk '
    NR==1 { for(i=1;i<=NF;i++) idx[$i]=i; next }
    { v = $(idx["LOG10P"]) + 0; if (NR==2 || v > max) { max=v; id=$(idx["ID"]); beta=$(idx["BETA"]); se=$(idx["SE"]); pval=$(idx["LOG10P"]) } }
    END { print id"\t"beta"\t"se"\t"pval }
  ' "$out_file")
  IFS=$'\t' read -r HIT_VARIANT HIT_BETA HIT_SE HIT_PVAL <<< "$result"
  echo "CANDIDATE VARIANT:$HIT_VARIANT"
  echo "Variant info from conditioned sumstats $(map_vals_to_string "$HIT_BETA" "$HIT_SE" "$HIT_PVAL")"
  if awk -v p="$HIT_PVAL" -v t="$threshold" 'BEGIN{exit !(p>t)}'; then
    NEXT_STEP=$((step+1))
  else
    NEXT_STEP=0
  fi
}

# ---- main: recursive conditional chain for one locus/region ----
run_locus() {
  local locus="$1" region_flag="$2"
  pretty_print "$locus $region_flag conditional chain." 50

  local log_file="${OUT}_${PHENO}_${locus}.log"
  echo "Logging to $log_file"
  local result_file="${OUT}_${PHENO}_${locus}.independent.snps"

  IFS=',' read -ra CONDITION_LIST <<< "$locus"
  local condition_variant="${CONDITION_LIST[-1]}"

  IFS=$'\t' read -r b0 s0 m0 <<< "$(get_sum_dict_data "$condition_variant")"
  {
    echo -e "VARIANT\tBETA\tSE\tMLOG10P\tBETA_cond\tSE_cond\tMLOG10P_cond\tVARIANT_cond"
    echo -e "${condition_variant}\t${b0}\t${s0}\t${m0}\tNA\tNA\tNA\tNA"
  } > "$result_file"
  echo "Variant info from original FG sumstats $(map_vals_to_string "$b0" "$s0" "$m0")"

  local step=1
  while (( step > 0 && step <= MAX_STEPS )); do
    regenie_run "$step" "$locus" "$region_flag" "$log_file"
    if [[ "$REGENIE_RET" -ne 0 ]]; then
      echo "RUN FAILED, check $log_file for errors." >&2
      return 1
    fi

    check_hit "$REGENIE_OUT_FILE" "$step" "$PVAL_THRESHOLD"
    step="$NEXT_STEP"
    if [[ "$step" -ne 0 ]]; then
      local cond_list_str; cond_list_str=$(IFS=,; echo "${CONDITION_LIST[*]}")
      IFS=$'\t' read -r cb cs cm <<< "$(get_sum_dict_data "$HIT_VARIANT")"
      echo -e "${HIT_VARIANT}\t${cb}\t${cs}\t${cm}\t${HIT_BETA}\t${HIT_SE}\t${HIT_PVAL}\t${cond_list_str}" >> "$result_file"
      condition_variant="$HIT_VARIANT"
      CONDITION_LIST+=("$HIT_VARIANT")
      echo "Variant info from original FG sumstats $(map_vals_to_string "$cb" "$cs" "$cm")"
      echo "Hit signficant, proceeding to condition..."
    else
      echo "Hit not significant. Ending loop."
    fi
  done
}

if [[ -n "$LOCUS_ARG" ]]; then
  check_region "$LOCUS_ARG" "$REGION_ARG"
  run_locus "$CHECKED_LOCUS" "$CHECKED_REGION_FLAG"
else
  while read -r loc reg; do
    [[ -n "$loc" ]] || continue
    check_region "$loc" "$reg"
    run_locus "$CHECKED_LOCUS" "$CHECKED_REGION_FLAG"
  done < "$LOCUS_LIST"
fi
