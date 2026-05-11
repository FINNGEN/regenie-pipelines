#!/usr/bin/awk -f
# Detects whether every non-NA row for the given phenotypes has the same sex value.
# Usage: zcat cov_pheno.gz | detect_single_sex.awk -v sexcol=SEX_IMPUTED -v phenocols=PHENO1,PHENO2
# Prints "1" (single sex) or "0" (mixed) to stdout. Errors on missing column.
BEGIN { is_single = 1 }
NR == 1 {
    for (i = 1; i <= NF; i++) h[$i] = i
    split(phenocols, ps, ",")
    prev = ""
    for (pi in ps) {
        if (!(ps[pi] in h)) {
            print "Given phenotype " ps[pi] " does not exist in phenofile" > "/dev/stderr"
            exit 1
        }
    }
    if (!(sexcol in h)) {
        print "Given sexcolumn:" sexcol " does not exist in phenofile" > "/dev/stderr"
        exit 1
    }
}
NR > 1 {
    for (pi in ps) {
        if ($h[ps[pi]] != "NA") {
            if (prev != "" && prev != $h[sexcol]) { is_single = 0; exit }
            prev = $h[sexcol]
        }
    }
}
END {
    print "is single ", is_single > "/dev/stderr"
    printf is_single
}
