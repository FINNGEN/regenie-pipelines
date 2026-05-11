#!/usr/bin/awk -f
# Post-processes a regenie .regenie.gz file that was produced from a pgen whose
# chromosome column still contains a "chr" prefix. Parses the SNPID (column ID,
# expected format like chrN_POS_REF_ALT) to recover the numeric chromosome and
# writes it back into column CHROM ($1). All other columns are preserved.
#
# Usage: zcat input.regenie.gz | strip_chr_from_chrom.awk | bgzip > output.regenie.gz
NR == 1 { print $0 }
NR > 1 {
    n_ = split($3, snpidlist, "_")
    gsub(/chr/, "", snpidlist[1])
    $1 = snpidlist[1]
    print $0
}
