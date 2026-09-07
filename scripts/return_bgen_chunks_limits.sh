#!/usr/bin/env bash
# Builds a manifest of bgen chunk genomic bounds (chrom, first/last variant position) by reading
# each chunk's local .bgi index via bgenix -list (index-only, no genotype data touched -- cheap even
# for hundreds of chunks). Runs in parallel across all local CPUs.
#
# Output paths are rewritten from the local scan path to their gs:// WDL-run equivalent (same
# CHROM.CHUNK naming, different root) -- so the manifest is built cheaply against an already-present
# local copy of the chunks, but the paths it emits are the real remote ones the WDL will localize.
#
# Usage: ./return_bgen_chunks_limits.sh [local_chunk_dir] [gcs_chunk_root] [out_tsv]
#   defaults: /home/pete/r14/bgen/release/data/chunk
#             gs://finngen-production-library-red/finngen_R14/bgen_1.0/data/chunk
#             bgen_chunks_manifest.tsv
#
# Output columns: PATH (gs://.../finngen_R14_CHROM.CHUNK.bgen)  CHROM (as reported by bgenix, e.g. chr1)  START_POS  END_POS
set -euo pipefail

LOCAL_DIR="${1:-/home/pete/r14/bgen/release/data/chunk}"
GCS_ROOT="${2:-gs://finngen-production-library-red/finngen_R14/bgen_1.0/data/chunk}"
OUT="${3:-bgen_chunks_manifest.tsv}"

chunk_bounds() {
  local f="$1" gcs_root="$2"
  local base chrom_chunk chrom chunk gcs_path
  base=$(basename "$f")
  # finngen_R14_<CHROM>.<CHUNK>.bgen -> CHROM, CHUNK
  chrom_chunk="${base#finngen_R14_}"
  chrom_chunk="${chrom_chunk%.bgen}"
  chrom="${chrom_chunk%%.*}"
  chunk="${chrom_chunk#*.}"
  gcs_path="${gcs_root}/finngen_R14_${chrom}.${chunk}.bgen"

  local vchrom start end
  read -r vchrom start end <<< "$(bgenix -g "$f" -list 2>/dev/null | awk '
    /^#/ {next}
    !h {h=1; next}
    {if (s=="") {s=$4; c=$3}; e=$4}
    END {print c, s, e}
  ')"
  printf "%s\t%s\t%s\t%s\n" "$gcs_path" "$vchrom" "$start" "$end"
}
export -f chunk_bounds

find "$LOCAL_DIR" -maxdepth 1 -name 'finngen_R14_*.*.bgen' \
  | parallel -j "$(nproc)" chunk_bounds {} "$GCS_ROOT" > "$OUT"

sort -t$'\t' -k2,2V -k3,3n -k1,1V -o "$OUT" "$OUT"

echo "Wrote $(wc -l < "$OUT") chunk records to $OUT" >&2
