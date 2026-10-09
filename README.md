# regenie-pipelines
WDL pipelines for running regenie

See regenie [documentation](https://rgcgithub.github.io/regenie/options/) and [paper](https://www.biorxiv.org/content/10.1101/2020.06.19.162354v2.full.pdf)

## Building docker image

The current version of regenie has the modifications created in the finngen repo, so it can be used to build the docker image. We want to build with Boost iostream for compression support, and with intel MKL as the linear algebra package.

Build image with [build_docker.sh](scripts/build_docker.sh). Give tag with version (e.g. FG_1.1) as parameter.
Check the beginning of scripts for modifying the default behavior with few variables (e.g. not building base regenie but basing off of already built one)
The parameter is just name tag to be added to the base regenie docker name. Change versioning as needed.


## Running GWAS

How to run regenie GWAS with Cromwell  
This in an example scenario creating new phenotypes with R7 data and running those

1. Create a covariate/phenotype file that contains your phenotypes. E.g. get `gs://r7_data/pheno/R7_COV_PHENO_V2.txt.gz`, add phenotypes to that (if a binary phenotype: cases 1, controls 0, everyone else NA - if a quantitative phenotype: inverse rank normalized values), and upload the new file to a bucket. Note that the phenotype file should be tab-separated with no spaces as both are treated as separator in regenie.
2. Create a text file with the names of your new phenotypes  
    2.1. If you have one or a few phenotypes, create a file with one phenotype per line, e.g.  
    my_phenos.txt
    ```
    PHENO1
    PHENO2
    ```
    and upload the file to a bucket.  
    2.2. If you have a larger number of phenotypes, create a tab-separated file where each line contains phenotypes that have < 5 % non-shared missingness. For example, if you have PHENO{1:5} that are female phenotypes (all males are NA and will be excluded as shared missingness) and they have < 5 % non-shared missingness among the females, and PHENO{6:7} that are male phenotypes with < 5 % non-shared missingness among the males, you can do:  
    my_phenos.txt
    ```
    PHENO1  PHENO2  PHENO3  PHENO4  PHENO5
    PHENO6  PHENO7
    ```
    and upload the file to a bucket. Phenotypes on each row will be analyzed together as regenie step 1 is faster that way. Non-shared missing phenotypes for phenotypes on each row will be mean-imputed for level 0 regression and this is why the phenotypes on each row should have low non-shared missingness. It's recommended to analyze at most 8 phenotypes together because the number of vCPUs used grows with the number of phenotypes and we've observed high preemption rates with VMs with more than 8 vCPUs.
3. Clone this repo: `git clone https://github.com/FINNGEN/regenie-pipelines`
4. Edit the input file `regenie-pipelines/wdl/gwas/regenie.json`:  
    5.1. Change `regenie.cov_pheno` to the file you created in the first step  
    5.2. Change `regenie.phenolist` to the file you created in the second step  
    5.3. `regenie.is_binary` should be `true` for binary phenotypes and `false` for quantitative phenotypes  
    5.4. `regenie.sub_step2.step2.test` can be `additive`, `recessive` or `dominant` depending on which analysis you want to run  
    5.5. Change `regenie.covariates` and `regenie.sub_step2.step2.options` as needed

    5.6.
    "regenie.auto_remove_sex_covar": true,
    "regenie.sex_col_name": "SEX_IMPUTED",
    "regenie.sub_step2.run_sex_specific": true,
    "regenie.sub_step2.step2.sex_specific_logpval": 6,

5. Cromwell requires subworkflows be zipped: `cd regenie-pipelines/wdl/gwas/ && zip regenie_sub regenie_step1.wdl regenie_sub.wdl`
6. Connect to Cromwell server  
    `gcloud compute ssh cromwell-fg-1 --project finngen-refinery-dev --zone europe-west1-b -- -fN -L localhost:5000:localhost:80`
7. Submit workflow  
    7.1. With `https://github.com/FINNGEN/CromwellInteract`  
    7.2. Or using the web interface  
        7.2.1 Go to `http://0.0.0.0:5000` with your browser  
        7.2.2 Click `/api/workflows/{version}`  
        7.2.3 Choose `regenie.wdl` as workflowSource  
        7.2.4 Choose the edited `regenie.json` as workflowInputs  
        7.2.5 Choose `regenie_sub.zip` as workflowDependencies  
        7.2.6 `Execute`
8. Use the given workflow id to look at the timing diagram or to get metadata  
`http://0.0.0.0:5000/api/workflows/v1/WORKFLOW_ID/timing`
`http://0.0.0.0:5000/api/workflows/v1/WORKFLOW_ID/metadata`
9. Logs and results go under  
`gs://fg-cromwell_fresh/regenie/WORKFLOW_ID`  
Summary stats and tabix indexes:  
`gs://fg-cromwell_fresh/regenie/WORKFLOW_ID/call-sub_step2/**/call-gather/**/*.gz*`
Plots:  
`gs://fg-cromwell_fresh/regenie/WORKFLOW_ID/call-sub_step2/**/*.png`  
Summary files with p < 1e-6 variants including annotation:  
`gs://fg-cromwell_fresh/regenie/WORKFLOW_ID/call-sub_step2/**/call-summary/**/*_summary.txt`


## FINNGEN CONDITIONAL ANALYSIS

This is a wrapper pipeline of [regenie](https://rgcgithub.github.io/regenie/) for conditional analysis. For each region it iteratively runs regenie step 2, each time conditioning on one more variant, until no significant hits are left. The wdl is meant for release purposes and can either discover all hits from a list of phenos based on the official FinnGen results, or run a user-supplied list of pheno/region/locus combinations.

The top-level [`regenie_conditional_analysis.wdl`](wdl/conditional-analysis/regenie_conditional_analysis.wdl) + [`regenie_conditional_analysis.json`](wdl/conditional-analysis/regenie_conditional_analysis.json) reflect the currently active version of the pipeline and are meant to be edited release over release. Because the wdl's input schema has changed between releases, each past release is frozen under its own `wdl/conditional-analysis/rXX/` folder, containing the exact wdl + json pair that release actually ran against — mirroring the per-release folder convention under `wdl/gwas/`, except here the wdl itself is versioned too since it isn't guaranteed stable across releases.

The sandbox (unmodifiable pipeline) version lives in the `sandbox-unmodifiable-pipelines` repo (`wdl/conditional/regenie_conditional_merged.sb.wdl`). It is an independently maintained copy of the same logic, so changes to the conditional logic here need to be ported there by hand.

### scripts/regenie_conditional.sh

This is the "engine" of the pipeline, that can also be used independently, so I will first explain its mechanism and inputs. The wdl does not call the script: its functions and driver logic are pasted verbatim into the `regenie_conditional` task's command block, with the WDL inputs assigned to the same shell variables the CLI parser would set. Any change to the script needs to be copied into the wdl (and vice versa).

These are the parameters:
```
Usage: regenie_conditional.sh --pheno P --out OUT --bgen B --sumstats S --null-file N
                               (--locus-region LOCUS REGION | --locus-list FILE) [options]
  --pval-threshold FLOAT     threshold limit (-log10(p)), or a raw p-value (<1) (default 7)
  --pheno-file FILE          pheno + covariate file
  --covariates LIST          comma-separated covariate list (default: full FinnGen list)
  --sample-file FILE         bgen sample file (auto-detected next to --bgen if omitted)
  --regenie-params STR       extra regenie params (default: " --bt --firth --approx --bsize 200 --ref-first")
  --force                    force re-run of already-completed steps
  --max-steps INT            default 10
  --chr-col/--pos-col/--ref-col/--alt-col/--mlogp-col/--beta-col/--sebeta-col STR
  --threads INT              default: nproc
```

They are all quite self explanatory. The null files are the `*loco.gz` outputs of regenie step 1. `--locus-region` and `--locus-list` are mutually exclusive and are meant for defining the regions of choice. The vanilla mode runs just one region/locus (in any order and in regenie format, e.g. `6:34869517-37869517 chr6_35376598_G_A`). The script will automatically recognize which is the locus and which the region. Else one can pass a file with a tsv separated list of regions/loci, one per line.

Each run will iteratively condition on more and more significant variants until no hits are found under a certain threshold (`--pval-threshold`, either a mlogp > 1 or a pval < 1, it gets converted to mlogp anyways). One can also cap the iterations at a certain number of steps (`--max-steps`) instead. The locus can also be a comma-separated list of variants, in which case all of them are conditioned on from the first step.

`--regenie-params` are the extra parameters to pass to regenie. `--ref-first` is required with FG data. For binary phenos Firth correction (`--firth --approx --pThresh ... --firth-se`) is recommended. The Firth null model is fit fresh at every step on purpose: reusing a null-firth file (either from step 1 or recycled across steps with `--write-null-firth`/`--use-null-firth`) was benchmarked to be up to ~7x *slower*, because the reused null has no estimate for the newest conditioning variant and regenie's solver falls into an expensive retry ladder. Don't reintroduce it without re-benchmarking.

The outputs will be in the `--out` directory (generated if missing). Along with a temporary folder that contains all the necessary files, the outputs are:
- prefix*_pheno_locus.log: the stdout/err of regenie is appended to this file so all logs are available
- prefix*_pheno_locus.independent.snps: contains the chain of results, with columns `VARIANT BETA SE MLOG10P BETA_cond SE_cond MLOG10P_cond VARIANT_cond`. The first row is the starting locus with its original sumstats values; each following row is a new independent hit with both its original and conditioned values, and `VARIANT_cond` lists the variants it was conditioned on.
- prefix*_pheno_locus_STEP.conditional: the regenie output of each of the [1..N] steps of the chain.

### WDL

Here I will explain the tasks and inputs of the [wdl](wdl/conditional-analysis/regenie_conditional_analysis.wdl)

#### Inputs

```
"regenie_conditional_analysis.docker": "eu.gcr.io/finngen-sandbox-v3-containers/regenie:4.1.2_cond_bgenix",
"regenie_conditional_analysis.test": false,
"regenie_conditional_analysis.pheno_region_input": "gs://.../phenos.txt",
"regenie_conditional_analysis.chroms": ["1", "2", ..., "22", "23"],
"regenie_conditional_analysis.release": "14",
"regenie_conditional_analysis.sumstats_root": "gs://r14-data/regenie/release/summary_stats/PHENO.gz",
"regenie_conditional_analysis.pheno_file": "gs://r14-data/pheno/R14_COV_PHENO_V0.FID.txt.gz",
"regenie_conditional_analysis.locus_mlogp_threshold": 7.3,
"regenie_conditional_analysis.conditioning_mlogp_threshold": 6,
"regenie_conditional_analysis.mlogp_col": "mlogp",
"regenie_conditional_analysis.chr_col": "#chrom",
"regenie_conditional_analysis.pos_col": "pos",
"regenie_conditional_analysis.ref_col": "ref",
"regenie_conditional_analysis.alt_col": "alt",
"regenie_conditional_analysis.chunk_manifest": "gs://finngen-production-library-green/wdl/conditional/bgen_chunks_manifest.tsv",
"regenie_conditional_analysis.covariates": ["SEX_IMPUTED", "AGE_AT_DEATH_OR_END_OF_FOLLOWUP", "PC1", ..., "PC10", "IS_FINNGEN2_CHIP", "BATCH_DS1_BOTNIA_Dgi_norm", ...],
```
`pheno_region_input` determines what is run, and its shape determines the mode (see `validate_regions` below). It's a tab separated file with no header and either:
- **1 column**: a list of phenos. Regions are discovered automatically from the finemap regions and sumstats (`extract_cond_regions` + `merge_regions`).
- **4 columns**: `pheno`, `chrom` (numeric, 23 for X), `region` (`chrom:start-end`), `locus` (variant ID(s) in `chrX_pos_ref_alt`-style format, comma-separated if more than one). Discovery is skipped entirely and these rows are run as they are. The same pheno can appear on several rows.

`release` sets the output prefix (`finngen_R<release>`).
`chroms` restricts region discovery to the given chromosomes. It's only applied in discovery mode, custom regions are run regardless.
`pheno_file` is the file that contains all pheno related data. It has to have the FID/IID columns and contain the covariates.
`sumstats_root` is the `PHENO`-templated path to the (tabixed!) sumstats; the column names are set by the `*_col` inputs.
`locus_mlogp_threshold` is the threshold for a region's top hit to be picked up as the starting point of a chain (discovery mode only).
`conditioning_mlogp_threshold` instead is the parameter used to stop the regenie chain conditional run.
`chunk_manifest` is the list of bgen chunks and their genomic bounds (see `attach_bgen_chunks` below).

`test`, when `true`, cuts the run down to a cheap smoke test: `pheno_region_input` is truncated to its first 10 rows in `validate_regions` (so 10 phenos in discovery mode, or 10 regions in custom mode), and in discovery mode `merge_regions` further keeps only 2 random regions per pheno. So a test run touches at most 10 phenos and 20 regions end to end.

There is no `is_binary` input anymore: binary/quantitative is detected automatically per pheno (see `check_is_binary`), so a batch can mix binary and quantitative phenos.

#### validate_regions
Checks the shape of `pheno_region_input`: all rows must have the same number of columns, and it must be either 1 or 4, otherwise the run fails. It returns which mode to run, the (possibly test-truncated) rows, and the deduplicated list of phenos (always the first column).

#### filter_covariates
This is a preprocessing task. It generates for each pheno the list of valid covariates to be passed to regenie. It checks that for each group of input phenos (in this case each pheno is its own group) there are at least N counts of non NA samples *and* non 0 covariates. The output of the task is a pheno --> covariates map object that is then passed to regenie later.
```
"regenie_conditional_analysis.filter_covariates.threshold_cov_count": 10,
```

#### check_is_binary
Classifies each pheno with a single pass over `pheno_file`: if all non-missing values (empty or `NA`) of a pheno's column are `0`/`1`, the pheno is binary, else it's quantitative. A pheno missing from the file fails the run. The result is a pheno --> is_binary map that `regenie_conditional` uses to pick between `regenie_params_binary` and `regenie_params_qt`. Note that this means binary phenos need to be coded 0/1 (a 1/2 coding would be treated as quantitative).

#### extract_cond_regions
Discovery mode only. This task returns the top hits for each pheno, above the `locus_mlogp_threshold`. The only required input are the finemap regions that are produced by our pipeline. For each region, the task filters the input sumstats (tabix file is to be expected!) to the region limits and returns the most significant variant, if it passes the threshold. If no region has a hit, the task will not fail: an empty file is generated regardless.

```
"regenie_conditional_analysis.extract_cond_regions.region_root": "gs://r14-data/finemap/release/beds/PHENO.bed",
"regenie_conditional_analysis.extract_cond_regions.add_hla": true,
```
`add_hla` appends the fixed HLA region (chr6:29,000,000-34,000,000) to every pheno's region bed before hit extraction, on top of whatever regions `region_root` already supplies.

The logic is the bash port in `scripts/filter_hits_regions.sh`, inlined in the task like `regenie_conditional.sh` is.

#### merge_regions
Discovery mode only. All regions from the previous task are merged into a single file with the same 4-column shape as the custom regions input, so the rest of the pipeline is the same for both modes.

#### attach_bgen_chunks
Instead of localizing the whole chromosome bgen for each region, only the bgen chunk(s) that overlap the region are localized. This task adds a 5th column to the regions file with the (comma separated) paths of the overlapping chunks, using `chunk_manifest`. The manifest is built once per release with [`scripts/return_bgen_chunks_limits.sh`](scripts/return_bgen_chunks_limits.sh), which reads each chunk's `.bgi` index and writes `path, chrom, start, end`. A region with no overlapping chunk fails the run.

#### regenie_conditional
This is the major task where the magic happens. The overlapping chunks are concatenated into a single local bgen (`cat-bgen`) and indexed, then the conditional chain is run on the region.
 ```
"regenie_conditional_analysis.regenie_conditional.null_root": "gs://r14-data/regenie/release/loco/R14_GRM_V0_LD_0.2.PHENO.loco.gz",
"regenie_conditional_analysis.regenie_conditional.beta": "beta",
"regenie_conditional_analysis.regenie_conditional.sebeta": "sebeta",
"regenie_conditional_analysis.regenie_conditional.max_steps": 10,
"regenie_conditional_analysis.regenie_conditional.regenie_params_binary": "--bt --firth --firth-se --approx --pThresh 0.01 --bsize 200 --ref-first",
"regenie_conditional_analysis.regenie_conditional.regenie_params_qt": "--qt --bsize 200 --ref-first",
"regenie_conditional_analysis.regenie_conditional.cpus": 4,
```
`null_root` are the step1 outputs.
`beta` and `sebeta` are the column names for the entries in the sumstat file.
`max_steps` controls the maximum length of the chain.
`cpus` is self explanatory.
The bgen chunks are expected to have a `.sample` file next to them (`<chunk>.bgen.sample`).

`regenie_params_binary`/`regenie_params_qt` are the full set of extra flags passed straight through to regenie, one of which is picked per pheno based on `check_is_binary` — there's no need to repeat `--bt`/`--qt` anywhere else. There is no null-firth file input: Firth is fit fresh at every step (see `regenie_conditional.sh` above for why).

#### merge_results
The chains are merged into one `independent_snps` file per pheno.

#### Outputs
- `all_chains`: the per-region chain files (`*.independent.snps`)
- `all_outputs`: the regenie output of every step (`*.conditional`)
- `pheno_chains`: the per-pheno merged chains

### PHEWEB IMPORT

The pheweb import munging is not part of the main wdl. There is a separate [`pheweb_import.wdl`](wdl/conditional-analysis/pheweb_import.wdl) that builds the sql import file and munges the regenie outputs, dealing with the files in chunks since the number of files in a release is large. The relevant inputs are the lists of paths for the conditional chains and the regenie outputs, plus the regions file, which are part of the release anyways.
