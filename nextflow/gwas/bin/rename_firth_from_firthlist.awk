#!/usr/bin/awk -f
# Emits `mv <orig> <renamed>` shell commands that rename regenie step1 firth files
# from `<prefix>_N.firth.gz` to `<prefix>.<pheno>.firth.gz`.
# Usage: rename_firth_from_firthlist.awk <prefix>_firth.list | bash
{ orig = $2; sub(/_[0-9]+.firth.gz/, "." $1 ".firth.gz", $2); print "mv " orig " " $2 }
