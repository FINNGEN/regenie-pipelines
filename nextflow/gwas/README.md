# FinnGen regenie GWAS — Nextflow DSL2 port

Nextflow DSL2 port of the FinnGen regenie GWAS WDL pipeline under
`wdl/gwas/`. Targets Google Cloud Batch. Supports three variants via a single
pipeline, selected by feature-flag params files:

| Variant                    | Params file                          | Genotype format | Summary/Coding | Sex-specific | Notes |
|----------------------------|--------------------------------------|-----------------|----------------|--------------|-------|
| Base release GWAS          | `params/r13_binary_release.json`     | bgen            | on             | on           | Mirrors `wdl/gwas/regenie.wdl` |
| Chip GWAS (bgen)           | `params/r13_chip_binary.json`        | bgen            | off            | off          | Mirrors `wdl/gwas/r13/chip/regenie.wdl` |
| Chip GWAS (pgen)           | `params/r13_chip_pgen.json`          | pgen            | off            | off          | Mirrors `wdl/gwas/r13/chip/pgen_pipeline/` — adds strip-chr post-processing, version-sort gather, disabled qqplot |

## Layout

```
nextflow/gwas/
├── main.nf                       # Workflow entrypoint
├── nextflow.config               # Profiles, resources, param defaults
├── modules/
│   ├── step1.nf                  # regenie step 1
│   ├── step2.nf                  # regenie step 2 (bgen + pgen)
│   ├── gather.nf                 # per-pheno concat + munge + qqplot
│   ├── summary.nf                # FinnGen annotation
│   └── coding_gather.nf          # Global coding-variant concat
├── bin/                          # Helper scripts (awk/python/bash)
└── params/                       # Per-variant JSON params files
```

## Running

```bash
# Stub run (smoke test — no real data needed):
nextflow run nextflow/gwas/main.nf -profile standard -stub \
    -params-file nextflow/gwas/params/r13_binary_release.json

# Real run on Google Cloud Batch:
nextflow run nextflow/gwas/main.nf -profile google_batch -resume \
    -params-file nextflow/gwas/params/r13_binary_release.json \
    -work-dir gs://<bucket>/nf-work
```

Each params file already sets `outdir` and `google_workdir` to reasonable
defaults — override on the command line if needed.

Use `-with-trace`, `-with-report`, `-with-timeline` for stats on the run.

## Feature flags

Set in the params JSON — all flags are documented in `nextflow.config`
and `main.nf`. The interesting ones:

| Flag                     | Meaning                                                  |
|--------------------------|----------------------------------------------------------|
| `input_format`           | `bgen` or `pgen` — selects the STEP2 genotype flags      |
| `step2_resource_profile` | `dynamic` (scale with pheno count) or `fixed`            |
| `step2_strip_chr`        | `true` to strip `chr` prefix from CHROM after STEP2      |
| `gather_sort_flag`       | sort key for the GATHER concat step                      |
| `skip_qqplot`            | `true` to skip the R qqplot in GATHER                    |
| `skip_summary`           | `true` to skip FinnGen annotation                        |
| `skip_coding_gather`     | `true` to skip the global coding-variant concat          |
| `run_sex_specific`       | `true` to run male/female stratified analysis in STEP2   |

## Notes on parity with the WDL

- The orphaned `gnomad_annotation` key from the original r13 JSONs is not
  consumed anywhere in the WDL (only `finngen_annotation` is). It has been
  dropped from the Nextflow params files.
- The `chip/pgen_pipeline` variant has a typo in its embedded annotation
  python (`EXOME_enrichment_nfsee`). Since chip variants set
  `skip_summary: true`, `annotate_summary.py` is never invoked there. The
  shared `bin/annotate_summary.py` uses the correct column names
  (`EXOME_enrichment_nfe`, `GENOME_enrichment_nfe`).
- The WDL chip-pgen `strip chr` post-processing only rewrites the *first*
  pheno's output in a chunk (`${phenolist[0]}`). The Nextflow port reproduces
  this behaviour faithfully for parity; if you want to strip every pheno's
  file, change the `strip_chr_post` block in `modules/step2.nf` to loop over
  all matching `.regenie.gz` files.
- Empty sex-specific shards are written with a valid header (see
  `bin/write_empty_sex_spec_header.sh`) so that GATHER can pick up the
  expected columns from the first shard, matching the WDL behaviour.
- STEP2 stages loco/nulls files with `stageAs: './*'` so regenie's
  `pred.list` / `firth.list` path references resolve against the task
  workdir — no need for the WDL's `mv $file .` loop.
- For the pgen variant, the params key is still named `bgenlist` for
  backwards compat with the original WDL JSONs. It should point at the
  pgen list file (e.g. `new_pgenlist`).
