#!/usr/bin/awk -f
# Convert a regenie binary output to the pheweb-import format (tab-separated TSV).
# Drops variants with LOG10P == NA (unsuccessful Firth/SPA).
# Usage: zcat regenie.gz | munge_regenie_binary.awk | bgzip > munged.gz
BEGIN {
    FS = " "
    OFS = "\t"
    split("CHROM GENPOS ALLELE0 ALLELE1 LOG10P BETA SE A1FREQ A1FREQ_CASES A1FREQ_CONTROLS", REQUIRED_FIELDS)
}
NR == 1 {
    for (i = 1; i <= NF; i++) h[$i] = i
    for (i in REQUIRED_FIELDS) {
        if (!(REQUIRED_FIELDS[i] in h)) {
            print REQUIRED_FIELDS[i] " expected in regenie header" >> "/dev/stderr"
            exit 1
        }
    }
    print "#chrom", "pos", "ref", "alt", "pval", "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls"
}
NR > 1 && $h["LOG10P"] != "NA" {
    print $h["CHROM"], $h["GENPOS"], $h["ALLELE0"], $h["ALLELE1"], 10^-$h["LOG10P"], $h["LOG10P"], $h["BETA"], $h["SE"], $h["A1FREQ"], $h["A1FREQ_CASES"], $h["A1FREQ_CONTROLS"]
}
