#!/usr/bin/awk -f
# Emits variant IDs with LOG10P above the given threshold from a regenie output.
# Usage: zcat <regenie.gz> | extract_sex_variants.awk -v threshold=1.623249 > variants.txt
NR == 1 {
    for (i = 1; i <= NF; i++) h[$i] = i
    if (!("ID" in h) || !("LOG10P" in h)) {
        print "ID and LOG10P columns expected in association file" > "/dev/stderr"
        exit 1
    }
}
NR > 1 && $h["LOG10P"] > threshold { print $h["ID"] }
