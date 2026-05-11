#!/usr/bin/awk -f
# Rewrites a firth.list so that entries `<prefix>_N.firth.gz` point to
# `<prefix>.<pheno>.firth.gz` (using col 1 as the pheno name).
# Usage: rewrite_firthlist.awk <prefix>_firth.list > <prefix>.<hash>.firth.list
{ sub(/_[0-9]+.firth.gz/, "." $1 ".firth.gz", $2); sub(/.*\//, "", $2) } 1
