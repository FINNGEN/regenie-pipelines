#!/usr/bin/awk -f
# Counts samples per sex matching a phenotype condition.
# For binary: counts rows where $pheno == 1.
# For quant : counts rows where $pheno != "NA".
# Samples are restricted to IDs in the given FID-IID keep file.
#
# Usage:
#   count_sex_cases.awk -v pheno=PHENO -v is_binary=1 males <(zcat cov_pheno.gz) | wc -l
#   count_sex_cases.awk -v pheno=PHENO -v is_binary=0 females <(zcat cov_pheno.gz) | wc -l
FNR == NR { keep[$1]; next }
FNR == 1  { for (i = 1; i <= NF; i++) h[$i] = i; next }
{
    if (is_binary == 1) {
        if ($h[pheno] == 1 && $1 in keep) print $1
    } else {
        if ($h[pheno] != "NA" && $1 in keep) print $1
    }
}
