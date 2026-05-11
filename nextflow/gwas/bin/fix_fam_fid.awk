#!/usr/bin/awk -f
# Sets FID column equal to IID column in a plink .fam file.
# Usage: fix_fam_fid.awk input.fam > patched.fam
{ $1 = $2 } 1
