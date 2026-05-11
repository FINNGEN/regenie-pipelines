#!/usr/bin/awk -f
# Filter binary covariates with fewer than `th` 0s or 1s among rows where at least
# one of the phenotypes is non-NA. Quantitative covariates and covariates not found
# in the header are passed through unchanged via a mask.
#
# Usage: zcat cov_pheno.gz | filter_covariates.awk -v covariates=COV1,COV2,PC{1:10} -v phenos=PHENO1,PHENO2 -v th=10
# Writes final comma-separated covariate list to stdout; logs per-covariate counts to stderr.
BEGIN { FS = "\t" }
NR == 1 {
    covlen = split(covariates, covars, ",")
    phlen = split(phenos, phenoarr, ",")
    for (i = 1; i <= NF; i++) {
        h[$i] = i
        mask[$i] = 0
    }
}
NR > 1 {
    # if any of the phenotypes is non-NA, take the row into account
    process_row = 0
    for (ph in phenoarr) {
        if ($h[phenoarr[ph]] != "NA") { process_row = 1 }
    }
    if (process_row == 1) {
        for (co in covars) {
            if ($h[covars[co]] == 0) {
                zerovals[covars[co]] += 1
            } else if ($h[covars[co]] == 1) {
                onevals[covars[co]] += 1
            } else if ($h[covars[co]] == "NA") {
                # no-op
                na = 0
            } else {
                # mask this covariate to be included, no matter the counts
                # includes both quantitative covars and covariates not found in header
                mask[covars[co]] = 1
            }
        }
    }
}
END {
    SEP = ""
    for (co in covars) {
        if ((zerovals[covars[co]] > th && onevals[covars[co]] > th) || mask[covars[co]] == 1) {
            printf("%s%s", SEP, covars[co])
            SEP = ","
        }
        printf "Covariate %s zero count: %d one count: %d mask: %d\n", covars[co], zerovals[covars[co]], onevals[covars[co]], mask[covars[co]] >> "/dev/stderr"
    }
}
