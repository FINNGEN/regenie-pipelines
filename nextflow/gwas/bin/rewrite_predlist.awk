#!/usr/bin/awk -f
# Rewrites a pred.list so that entries `<prefix>_N.loco.gz` point to
# `<prefix>.<pheno>.loco.gz` (using col 1 as the pheno name).
# Usage: rewrite_predlist.awk <prefix>_pred.list > <prefix>.<hash>.pred.list
{ sub(/_[0-9]+.loco.gz/, "." $1 ".loco.gz", $2); sub(/.*\//, "", $2) } 1
