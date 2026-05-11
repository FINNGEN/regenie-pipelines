#!/usr/bin/awk -f
# Emits `mv <orig> <renamed>` shell commands that rename regenie step1 loco files
# from `<prefix>_N.loco.gz` to `<prefix>.<pheno>.loco.gz`, where the pheno name is
# taken from column 1 of the pred.list file.
# Usage: rename_loco_from_predlist.awk <prefix>_pred.list | bash
{ orig = $2; sub(/_[0-9]+.loco.gz/, "." $1 ".loco.gz", $2); print "mv " orig " " $2 }
