#!/usr/bin/env python3
"""Merge male- and female-only regenie results with the base (both-sex) results,
and compute a per-variant male/female effect-size difference z-test.

Usage: merge_sex_specific.py <pheno> <prefix> <base_regenie_gz>

Inputs:
    pheno              Phenotype name used in the regenie output filenames.
    prefix             STEP2 filename prefix (same as in the WDL).
    base_regenie_gz    Full path to the base (both-sex) regenie output for this pheno.

Expects the following files to exist in the working directory:
    <prefix>.sex_spec.males_<pheno>.gz
    <prefix>.sex_spec.females_<pheno>.gz

Writes:
    <prefix>.<pheno>.sex_spec.gz
"""
import math
import sys

import gzip
import numpy as np
import pandas as pd
from scipy.stats import norm


def main(pheno: str, prefix: str, base: str) -> None:
    male = f"{prefix}.sex_spec.males_{pheno}.gz"
    female = f"{prefix}.sex_spec.females_{pheno}.gz"

    # NOTE: See write_empty_sex_spec_header.sh — we echo the expected header for empty
    # shards so that the gather step can take the header from the first shard. If you
    # change the output columns here, update that helper to match.
    sex_cols = ["ID", "N", "BETA", "SE", "LOG10P"]
    binarycols = ["A1FREQ_CONTROLS", "A1FREQ_CASES"]

    basic = pd.read_csv(gzip.open(base), sep=" ")
    cols = list(basic.columns)

    malestat = pd.read_csv(gzip.open(male), sep=" ")
    femalestat = pd.read_csv(gzip.open(female), sep=" ")

    if all(c in cols for c in binarycols):
        sex_cols.extend(binarycols)

    malestat = malestat[sex_cols]
    femalestat = femalestat[sex_cols]

    combs = (
        basic.merge(malestat.add_prefix("males_"), left_on="ID", right_on="males_ID")
        .merge(femalestat.add_prefix("females_"), left_on="ID", right_on="females_ID")
    )

    combs["diff_beta"] = combs["males_BETA"] - combs["females_BETA"]

    if len(combs.index) > 0:
        combs["p_diff"] = -(
            (
                norm.logsf(
                    abs(combs["diff_beta"]) / np.sqrt(combs["males_SE"] ** 2 + combs["females_SE"] ** 2)
                )
                + math.log(2)
            )
            / math.log(10)
        )
    else:
        combs["p_diff"] = 0

    combs.to_csv(f"{prefix}.{pheno}.sex_spec.gz", compression="gzip", sep=" ", index=False)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: merge_sex_specific.py <pheno> <prefix> <base_regenie_gz>", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
