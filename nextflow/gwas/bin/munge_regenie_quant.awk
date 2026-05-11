#!/usr/bin/awk -f
# Convert a regenie quantitative output to the pheweb-import format (tab-separated TSV).
# Usage: zcat regenie.gz | munge_regenie_quant.awk | bgzip > munged.gz
BEGIN {
    FS = " "
    OFS = "\t"
    split("CHROM GENPOS ALLELE0 ALLELE1 LOG10P BETA SE A1FREQ", REQUIRED_FIELDS)
}
NR == 1 {
    for (i = 1; i <= NF; i++) h[$i] = i
    for (i in REQUIRED_FIELDS) {
        if (!(REQUIRED_FIELDS[i] in h)) {
            print REQUIRED_FIELDS[i] " expected in regenie header" >> "/dev/stderr"
            exit 1
        }
    }
    print "#chrom", "pos", "ref", "alt", "pval", "mlogp", "beta", "sebeta", "af_alt"
}
NR > 1 {
    print $h["CHROM"], $h["GENPOS"], $h["ALLELE0"], $h["ALLELE1"], 10^-$h["LOG10P"], $h["LOG10P"], $h["BETA"], $h["SE"], $h["A1FREQ"]
}
