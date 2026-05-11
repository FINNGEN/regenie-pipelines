#!/usr/bin/awk -f
# Prints "FID IID" pairs for rows whose sex column equals the requested value.
# Usage: zcat cov_pheno.gz | extract_sex_lists.awk -v sexcol=SEX_IMPUTED -v sexval=0 > males
#        zcat cov_pheno.gz | extract_sex_lists.awk -v sexcol=SEX_IMPUTED -v sexval=1 > females
NR == 1 {
    for (i = 1; i <= NF; i++) h[$i] = i
    if (!(sexcol in h)) {
        print "Given sex column " sexcol " not found in phenotype file" > "/dev/stderr"
        exit 1
    }
}
NR > 1 && $h[sexcol] == sexval { print $1, $1 }
