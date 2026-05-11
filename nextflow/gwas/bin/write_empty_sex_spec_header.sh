#!/usr/bin/env bash
# Writes a header-only bgzipped sex_spec file. GATHER takes the header from the
# first shard, so empty shards must still contain a valid header line.
#
# Usage: write_empty_sex_spec_header.sh <is_binary:true|false> <out_file>
set -euo pipefail

is_binary=${1:?"is_binary arg required"}
out_file=${2:?"out_file arg required"}

if [[ "$is_binary" == "true" ]]; then
    echo -n "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ A1FREQ_CASES A1FREQ_CONTROLS INFO N TEST BETA SE CHISQ LOG10P"\
" EXTRA males_ID males_A1FREQ_CASES males_A1FREQ_CONTROLS males_N males_BETA males_SE males_LOG10P"\
" females_ID females_A1FREQ_CASES females_A1FREQ_CONTROLS females_N females_BETA females_SE females_LOG10P"\
" diff_beta p_diff" | bgzip > "$out_file"
else
    echo -n "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P"\
" EXTRA males_ID males_N males_BETA males_SE males_LOG10P"\
" females_ID females_N females_BETA females_SE females_LOG10P"\
" diff_beta p_diff" | bgzip > "$out_file"
fi
