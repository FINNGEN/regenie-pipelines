version 1.0

# Same workflow as conditional_merged.wdl, except regenie_conditional runs the bash port
# (scripts/regenie_conditional.sh) inlined directly in the command block instead of shelling out to
# regenie_conditional.py. WDL inputs are assigned straight into the shell variables the script's own
# functions expect (PHENO, OUT, BGEN, ...), then the script's functions/driver logic are pasted in as-is.
# Keep in sync with scripts/regenie_conditional.sh.

workflow conditional_analysis_bash {

  input {
    String docker
    # 1 column (pheno) -> region discovery runs via extract_cond_regions/merge_regions.
    # 4 columns (pheno,chrom,chrom:start-end,locus) -> treated as pre-supplied custom regions, discovery is skipped.
    File pheno_region_input
    String release
    Boolean test
    String mlogp_col
    String chr_col
    String pos_col
    String ref_col
    String alt_col
    Float conditioning_mlogp_threshold
    Float locus_mlogp_threshold
    Array[String] chroms
    Array[String] covariates
    File pheno_file
    String sumstats_root
    File chunk_manifest
  }

  String prefix = "finngen_R" + release

  # figure out which shape pheno_region_input is, and apply test-mode truncation (head -n 10) once,
  # uniformly, here -- so test mode runs either 10 phenos (discovery) or 10 regions (custom), not both.
  call validate_regions {input: pheno_region_input=pheno_region_input, test=test}
  Boolean is_custom_regions = validate_regions.is_custom_regions
  Array[String] pheno_data = validate_regions.pheno_data

  # returns covariate string for each pheno
  call filter_covariates {input: pheno_file=pheno_file,pheno_list=write_lines(pheno_data),covariates=covariates}
  # returns is_binary (0/1) for each pheno, one awk pass over the whole batch instead of per-pheno
  call check_is_binary {input: pheno_file=pheno_file,pheno_list=write_lines(pheno_data)}
  Map[String,String] is_binary_map = read_map(check_is_binary.is_binary_tsv)
  Map[String,String] cov_map = read_map(filter_covariates.cov_pheno_map)

  # region discovery -- skipped entirely when pheno_region_input was already a custom regions file
  if (!is_custom_regions) {
    # go through phenos
    scatter (p in pheno_data) {
      #get hits under pval threshold
      call extract_cond_regions {
        input: pheno=p,mlogp_threshold=locus_mlogp_threshold,mlogp_col=mlogp_col,chr_col=chr_col,pos_col=pos_col,ref_col=ref_col,alt_col=alt_col,chroms=chroms,sumstats_root=sumstats_root
      }
    }

    #SINGLE JOB PER LOCUS VERSION
    call merge_regions {input: hits=extract_cond_regions.gw_sig_res, test=test}
  }

  # either the user-supplied regions or the discovered ones -- same 4-column (pheno,chrom,region,locus)
  # shape either way, so attach_bgen_chunks runs once on the unified file regardless of branch.
  File region_rows_4col = if is_custom_regions then validate_regions.rows else select_first([merge_regions.regions])
  call attach_bgen_chunks {input: regions=region_rows_4col, chunk_manifest=chunk_manifest}
  Array[Array[String]] all_regions = read_tsv(attach_bgen_chunks.regions_with_chunks)

  # loop over all regions running one variant per shard
  scatter (region in all_regions) {
    String pheno = region[0]
    String chrom = region[1]
    String region_limits = region[2]
    String locus = region[3]
    Boolean pheno_is_binary = is_binary_map[pheno] == "1"

    # WDL 1.0 has no split() builtin: replace commas with newlines, then round-trip through
    # write_lines/read_lines to get a real Array[String] -- which then coerces elementwise to
    # Array[File], same as this WDL's existing String->File tricks (e.g. sub(sumstats_root,...)).
    Array[String] bgen_chunk_paths = read_lines(write_lines([sub(region[4], ",", "\n")]))
    Array[File] bgen_chunks = bgen_chunk_paths

    call regenie_conditional_bash {
      input: docker=docker,prefix=prefix,locus=locus,region=region_limits,pheno=pheno,chrom=chrom,covariates=cov_map[pheno],mlogp_col=mlogp_col,chr_col=chr_col,pos_col=pos_col,ref_col=ref_col,alt_col=alt_col,pval_threshold=conditioning_mlogp_threshold,sumstats_root=sumstats_root,pheno_file=pheno_file,is_binary=pheno_is_binary,bgen_chunks=bgen_chunks
    }
  }

  Array[File] results = flatten(regenie_conditional_bash.conditional_chains)
  call merge_results {input: prefix=prefix,result_list=results,phenos_list=write_lines(pheno_data)}

  output {
    Array[File] all_chains = flatten(regenie_conditional_bash.conditional_chains)
    Array[File] all_outputs = flatten(regenie_conditional_bash.regenie_output)
    Array[File] pheno_chains = merge_results.pheno_independent_snps
  }
}

task validate_regions {

  input {
    File pheno_region_input
    Boolean test
  }

  Int disk_size = ceil(size(pheno_region_input,'GB')) + 2

  command <<<
    set -eu
    cp ~{pheno_region_input} raw_input.tsv

    # single truncation point for test mode: 10 phenos-to-discover, or 10 pre-supplied regions -- never both
    if ~{true="true" false="false" test}; then
      head -n 10 raw_input.tsv > rows.tsv
    else
      cp raw_input.tsv rows.tsv
    fi

    ncols=$(awk -F'\t' '{print NF}' rows.tsv | sort -u)
    if [[ $(echo "$ncols" | wc -l) -ne 1 ]]; then
      echo "ERROR: pheno_region_input has inconsistent column counts across rows: $(echo "$ncols" | tr '\n' ' ')" >&2
      exit 1
    fi
    case "$ncols" in
      1) echo false > is_custom.txt ;;
      4) echo true > is_custom.txt ;;
      *) echo "ERROR: pheno_region_input must have 1 column (phenos to discover regions for) or 4 columns (pheno,chrom,region,locus custom regions), got $ncols" >&2; exit 1 ;;
    esac

    # pheno is always column 1 regardless of shape; dedupe once here since a pheno can appear on multiple
    # rows in custom-regions mode (e.g. several loci for the same pheno) -- downstream consumers don't need to.
    cut -f1 rows.tsv | sort -u > pheno_data.txt
  >>>

  output {
    Boolean is_custom_regions = read_boolean("is_custom.txt")
    File rows = "rows.tsv"
    Array[String] pheno_data = read_lines("pheno_data.txt")
  }

  runtime {
    disks: "local-disk ${disk_size} HDD"
  }
}

task merge_results {

  input {
    Array[File] result_list
    File phenos_list
    String prefix
  }

  command <<<
    cat ~{write_lines(result_list)} > results.txt
    # add prefix and suffix to make sure matching is not based on partially similar pheno names (ICD roots)
    cat ~{phenos_list} | awk '{print "'~{prefix}'_"$0"_chr"}' > phenos.txt

    # match pheno names and files and extract phenos back by removing pre/suffix
    grep -of phenos.txt results.txt | sed 's/_chr//g' > match_phenos.txt
    # same grepping but only keeping file path
    grep -f phenos.txt results.txt > match_paths.txt

    #now paste together the two files to have a tab separated pheno and filepath
    paste match_phenos.txt match_paths.txt > tmp.txt && head tmp.txt

    # write header in each pheno file
    head -n1 results.txt | xargs head -n1 > header.txt
    # loop through phenos and write header line to each pheno output
    while read f; do cp header.txt $f.independent_snps.txt; done < <(cut -f1 tmp.txt)
    # append to pheno output the content of each matching regenie output
    while read line; do arr=($line) && echo ${arr[0]} && cat ${arr[1]} | sed -e 1d >> ${arr[0]}.independent_snps.txt; done < tmp.txt
  >>>

  output {
    Array[File] pheno_independent_snps = glob("${prefix}*.txt")
  }
}

task attach_bgen_chunks {
  # Adds a 5th column (comma-joined gs:// chunk paths overlapping the region) to the 4-column
  # (pheno, chrom, region, locus) regions file, using a chunk-bounds manifest built by
  # scripts/return_bgen_chunks_limits.sh. Runs once, on the unified regions file, regardless of
  # whether it came from region discovery or custom user-supplied regions -- so both paths get
  # chunk-aware bgen loading identically.
  #
  # Deliberately outputs ONE combined file rather than one-file-per-region: with ~16k regions in
  # a full production run, glob-ing thousands of tiny per-region outputs would mean Cromwell
  # delocalizing thousands of files sequentially (observed ~1-2s fixed overhead per file on this
  # backend even for trivial files) -- hours of pure bookkeeping before any real work starts. A
  # single combined TSV is one delocalization; each scatter shard splits its own row's chunk-path
  # column into an Array[File] itself (see workflow body).

  input {
    File regions
    File chunk_manifest
  }

  String outfile = "regions_with_chunks.tsv"

  command <<<
    set -euo pipefail
    awk -F'\t' '
      FNR==NR {
        n[$2]++
        path[$2,n[$2]] = $1
        start[$2,n[$2]] = $3
        end[$2,n[$2]] = $4
        next
      }
      {
        split($3, rc, ":"); split(rc[2], se, "-")
        rstart = se[1]; rend = se[2]
        # matches the same 23<->X mapping used for variant IDs elsewhere: the manifest'\''s chrom
        # label comes straight from bgenix, which reports the X bgen'\''s chromosome as "X", not "23"
        key = ($2 == "23") ? "chrX" : "chr" $2
        chunks = ""
        for (i=1; i<=n[key]; i++) {
          if (end[key,i] >= rstart && start[key,i] <= rend) {
            chunks = (chunks == "") ? path[key,i] : chunks "," path[key,i]
          }
        }
        # never leave chunks empty -- a trailing empty TSV field gets silently dropped by
        # Cromwell'\''s read_tsv (confirmed: turns a 5-column row into 4, breaking region[4] downstream
        # with an out-of-bounds error instead of a clear failure). A region with zero overlapping
        # chunks means a real manifest gap; fail loudly here instead of producing a bad row.
        if (chunks == "") {
          print "ERROR: no bgen chunks overlap region " $3 " (chrom key " key ")" > "/dev/stderr"
          exit 1
        }
        print $1 "\t" $2 "\t" $3 "\t" $4 "\t" chunks
      }
    ' ~{chunk_manifest} ~{regions} > "~{outfile}"
  >>>

  output {
    File regions_with_chunks = outfile
  }
}

task regenie_conditional_bash {

  input {
    # GENERAL PARAMS
    String? regenie_docker
    String docker
    String prefix
    # hit info
    String locus
    String region
    String pheno
    String chrom
    # files to localize
    File pheno_file
    Array[File] bgen_chunks
    String null_root
    String sumstats_root
    # column names and stuff
    String chr_col
    String pos_col
    String ref_col
    String alt_col
    String mlogp_col
    String beta
    String sebeta
    # Script parameters/options
    Float pval_threshold
    Int max_steps
    String covariates
    Boolean is_binary
    String regenie_params_binary
    String regenie_params_qt
    Int cpus
  }

  # localize all files based on roots and pheno/chrom info
  File sumstats = sub(sumstats_root,"PHENO",pheno)
  File sum_tabix = sumstats + ".tbi"
  File null = sub(null_root,"PHENO",pheno)
  # all chunks for a chromosome share the same sample set/order -- any one chunk's .sample file works
  File bgen_sample = bgen_chunks[0] + ".sample"

  # is_binary picks which json-supplied param string applies
  String selected_params = if is_binary then regenie_params_binary else regenie_params_qt

  # runtime params based on file sizes
  Int disk_size = 120
  String final_docker = if defined(regenie_docker) then regenie_docker else docker

  command <<<
    set -euo pipefail
    echo ~{pheno} ~{chrom} ~{cpus}
    tabix -h ~{sumstats} ~{region} > region_sumstats.txt

    # merge only the chunk(s) overlapping this region into one local bgen, instead of localizing
    # the whole chromosome -- chunks need no pre-existing .bgi (cat-bgen concatenates raw variant
    # blocks directly), so only the merged result gets indexed.
    cat-bgen -og ./region.bgen -clobber -g ~{sep=' ' bgen_chunks}
    bgenix -g ./region.bgen -index -clobber

    # ---- WDL inputs -> the shell variables scripts/regenie_conditional.sh's own CLI parser would set ----
    PHENO="~{pheno}"
    OUT="./~{prefix}"
    PHENO_FILE="~{pheno_file}"
    BGEN="./region.bgen"
    SAMPLE_FILE="~{bgen_sample}"
    SUMSTATS="region_sumstats.txt"
    NULL_FILE="~{null}"
    COVARIATES="~{covariates}"
    REGENIE_PARAMS="~{selected_params}"
    REGENIE_CMD="regenie"
    FORCE=0
    MAX_STEPS=~{max_steps}
    CHR_COL="~{chr_col}"
    POS_COL="~{pos_col}"
    REF_COL="~{ref_col}"
    ALT_COL="~{alt_col}"
    MLOGP_COL="~{mlogp_col}"
    BETA_COL="~{beta}"
    SEBETA_COL="~{sebeta}"
    THREADS=~{cpus}
    PVAL_THRESHOLD=~{pval_threshold}
    LOCUS_ARG="~{locus}"
    REGION_ARG="~{region}"

    # ---- everything below is scripts/regenie_conditional.sh's functions + driver logic, pasted verbatim
    #      (CLI parsing/validation dropped -- the variables above already replace it). Keep in sync. ----

    pretty_print() {
      local string="$1" l="${2:-30}"
      local half=$(( l - ${#string}/2 )); (( half < 0 )) && half=0
      local dashes; dashes=$(printf '%*s' "$half" '' | tr ' ' '-')
      echo "${dashes}> ${string} <${dashes}"
    }

    map_vals_to_string() {
      awk -v b="$1" -v s="$2" -v m="$3" 'BEGIN{printf "{'"'"'beta'"'"': %.2f, '"'"'sebeta'"'"': %.2f, '"'"'mlogp'"'"': %.2f}", b, s, m}'
    }

    # region[0] can be "locus" or "region" depending on argument order -- whichever looks like CHR:START-END is the region.
    check_region() {
      local locus="$1" region="$2"
      if [[ "$locus" == *:* && "$locus" == *-* ]]; then
        local tmp="$locus"; locus="$region"; region="$tmp"
      fi
      CHECKED_LOCUS="$locus"
      CHECKED_REGION_FLAG=" --range $region "
    }

    OUT_DIR=$(dirname "$OUT")
    BASENAME=$(basename "$OUT")
    TMP_DIR="$OUT_DIR/tmp"
    mkdir -p "$OUT_DIR" "$TMP_DIR"

    # automatically look for sample file if not explicitly passed
    if [[ -z "$SAMPLE_FILE" ]]; then
      for candidate in "${BGEN%.bgen}.bgen.sample" "${BGEN%.bgen}.sample"; do
        [[ -f "$candidate" ]] && { SAMPLE_FILE="$candidate"; echo "WARNING: Using $SAMPLE_FILE as sample file." >&2; break; }
      done
    fi

    # figure out whether threshold is a p-value or -log10(p)
    PVAL_THRESHOLD=$(awk -v v="$PVAL_THRESHOLD" 'BEGIN{printf "%.10g", (v<1) ? -log(v)/log(10) : v}')
    pretty_print "MLOGP THRESHOLD: $PVAL_THRESHOLD"

    # ---- sumstats lookup cache: snpid -> beta, sebeta, mlogp (from the ORIGINAL, unconditioned sumstats) ----
    PVAL_CACHE="$TMP_DIR/${PHENO}_pvals.tsv"
    if [[ ! -f "$PVAL_CACHE" || "$FORCE" -eq 1 ]]; then
      echo "reading original pvals.." >&2
      zcat -f "$SUMSTATS" | awk -F'\t' -v chr_col="$CHR_COL" -v pos_col="$POS_COL" -v ref_col="$REF_COL" -v alt_col="$ALT_COL" -v mlogp_col="$MLOGP_COL" -v beta_col="$BETA_COL" -v sebeta_col="$SEBETA_COL" '
        NR==1 {
          for (i=1;i<=NF;i++) idx[$i]=i
          c=idx[chr_col]; p=idx[pos_col]; r=idx[ref_col]; a=idx[alt_col]; m=idx[mlogp_col]; b=idx[beta_col]; s=idx[sebeta_col]
          next
        }
        { print "chr" $c "_" $p "_" $r "_" $a "\t" $b "\t" $s "\t" $m }
      ' > "$PVAL_CACHE"
      echo "done." >&2
    fi

    # get_sum_dict_data equivalent: prints "beta\tsebeta\tmlogp" for a variant, or "0\t0\t0" if not found
    # (matches the python script's defaultdict(lambda: defaultdict(lambda: "0")) fallback)
    get_sum_dict_data() {
      local variant="$1"
      local line
      line=$(awk -F'\t' -v v="$variant" '$1==v{print $2"\t"$3"\t"$4; exit}' "$PVAL_CACHE")
      [[ -n "$line" ]] && echo "$line" || echo -e "0\t0\t0"
    }

    # ---- filter_pheno: extract FID/IID + pheno + covariate columns into a small tmp pheno file ----
    FILTERED_PHENO_FILE="$TMP_DIR/${PHENO}_pheno.tmp"
    if [[ ! -f "$FILTERED_PHENO_FILE" || "$FORCE" -eq 1 ]]; then
      echo "Creating new pheno file..." >&2
      header=$(zcat -f "$PHENO_FILE" | { head -1; cat >/dev/null; })
      cols=$(awk -F'\t' -v want="${PHENO},${COVARIATES}" -v h="$header" '
        BEGIN {
          n = split(want, w, ",")
          split(h, f, "\t")
          for (i in f) idx[f[i]] = i
          out = "1,2"
          for (j=1; j<=n; j++) {
            if (!(w[j] in idx)) { print "column " w[j] " not found in pheno file" > "/dev/stderr"; exit 1 }
            out = out "," idx[w[j]]
          }
          print out
        }')
      zcat -f "$PHENO_FILE" | cut -f"$cols" > "$FILTERED_PHENO_FILE"
      echo "done." >&2
    fi

    # ---- single regenie conditional run; sets REGENIE_OUT_FILE and REGENIE_RET ----
    regenie_run() {
      local step="$1" locus="$2" region_flag="$3" log_file="$4"

      local regenie_file="$TMP_DIR/${BASENAME}_${PHENO}.regenie"
      local out_file="$OUT_DIR/${BASENAME}_${PHENO}_${locus}_${step}.conditional"
      pretty_print "VARIANT:${CONDITION_LIST[-1]}"
      echo "generating $out_file ..." >&2

      if [[ ! -f "$out_file" || "$FORCE" -eq 1 ]]; then
        FORCE=1

        local pred_file="$TMP_DIR/${BASENAME}_${PHENO}.pred"
        printf "%s\t%s\n" "$PHENO" "$NULL_FILE" > "$pred_file"

        local tmp_variant="$TMP_DIR/${BASENAME}_variant.tmp"
        printf "%s\n" "${CONDITION_LIST[@]}" > "$tmp_variant"

        local sample_cmd=""
        [[ -n "$SAMPLE_FILE" && -f "$SAMPLE_FILE" ]] && sample_cmd=" --sample $SAMPLE_FILE"

        local cmd="$REGENIE_CMD --step 2 $REGENIE_PARAMS --bgen $BGEN${sample_cmd} --out $TMP_DIR/$BASENAME --pred $pred_file --phenoFile $FILTERED_PHENO_FILE --phenoCol $PHENO --condition-list $tmp_variant $region_flag --covarFile $FILTERED_PHENO_FILE --covarColList $COVARIATES --threads $THREADS"
        echo "$cmd" >&2

        local start=$SECONDS ret=0
        bash -c "$cmd" >> "$log_file" || ret=$?
        echo "Script ran in $((SECONDS-start)) seconds with $THREADS cpus." >&2

        [[ "$ret" -eq 0 ]] && mv "$regenie_file" "$out_file"
        REGENIE_OUT_FILE="$out_file"
        REGENIE_RET="$ret"
      else
        echo "file already exists" >&2
        REGENIE_OUT_FILE="$out_file"
        REGENIE_RET=0
      fi
    }

    # ---- check_hit: find the top (max LOG10P) row in a regenie output file; sets NEXT_STEP and HIT_* ----
    check_hit() {
      local out_file="$1" step="$2" threshold="$3"
      local result
      result=$(awk '
        NR==1 { for(i=1;i<=NF;i++) idx[$i]=i; next }
        { v = $(idx["LOG10P"]) + 0; if (NR==2 || v > max) { max=v; id=$(idx["ID"]); beta=$(idx["BETA"]); se=$(idx["SE"]); pval=$(idx["LOG10P"]) } }
        END { print id"\t"beta"\t"se"\t"pval }
      ' "$out_file")
      IFS=$'\t' read -r HIT_VARIANT HIT_BETA HIT_SE HIT_PVAL <<< "$result"
      echo "CANDIDATE VARIANT:$HIT_VARIANT"
      echo "Variant info from conditioned sumstats $(map_vals_to_string "$HIT_BETA" "$HIT_SE" "$HIT_PVAL")"
      if awk -v p="$HIT_PVAL" -v t="$threshold" 'BEGIN{exit !(p>t)}'; then
        NEXT_STEP=$((step+1))
      else
        NEXT_STEP=0
      fi
    }

    # ---- main: recursive conditional chain for one locus/region ----
    run_locus() {
      local locus="$1" region_flag="$2"
      pretty_print "$locus $region_flag conditional chain." 50

      local log_file="${OUT}_${PHENO}_${locus}.log"
      echo "Logging to $log_file"
      local result_file="${OUT}_${PHENO}_${locus}.independent.snps"

      IFS=',' read -ra CONDITION_LIST <<< "$locus"
      local condition_variant="${CONDITION_LIST[-1]}"

      IFS=$'\t' read -r b0 s0 m0 <<< "$(get_sum_dict_data "$condition_variant")"
      {
        echo -e "VARIANT\tBETA\tSE\tMLOG10P\tBETA_cond\tSE_cond\tMLOG10P_cond\tVARIANT_cond"
        echo -e "${condition_variant}\t${b0}\t${s0}\t${m0}\tNA\tNA\tNA\tNA"
      } > "$result_file"
      echo "Variant info from original FG sumstats $(map_vals_to_string "$b0" "$s0" "$m0")"

      local step=1
      while (( step > 0 && step <= MAX_STEPS )); do
        regenie_run "$step" "$locus" "$region_flag" "$log_file"
        if [[ "$REGENIE_RET" -ne 0 ]]; then
          echo "RUN FAILED, check $log_file for errors." >&2
          return 1
        fi

        check_hit "$REGENIE_OUT_FILE" "$step" "$PVAL_THRESHOLD"
        step="$NEXT_STEP"
        if [[ "$step" -ne 0 ]]; then
          local cond_list_str; cond_list_str=$(IFS=,; echo "${CONDITION_LIST[*]}")
          IFS=$'\t' read -r cb cs cm <<< "$(get_sum_dict_data "$HIT_VARIANT")"
          echo -e "${HIT_VARIANT}\t${cb}\t${cs}\t${cm}\t${HIT_BETA}\t${HIT_SE}\t${HIT_PVAL}\t${cond_list_str}" >> "$result_file"
          condition_variant="$HIT_VARIANT"
          CONDITION_LIST+=("$HIT_VARIANT")
          echo "Variant info from original FG sumstats $(map_vals_to_string "$cb" "$cs" "$cm")"
          echo "Hit signficant, proceeding to condition..."
        else
          echo "Hit not significant. Ending loop."
        fi
      done
    }

    check_region "$LOCUS_ARG" "$REGION_ARG"
    run_locus "$CHECKED_LOCUS" "$CHECKED_REGION_FLAG"
  >>>

  output {
    Array[File] conditional_chains = glob("./${prefix}*.snps")
    Array[File] logs = glob("./${prefix}*.log")
    Array[File] regenie_output = glob("./${prefix}*.conditional")
  }

  runtime {
    cpu: "~{cpus}"
    docker: "${final_docker}"
    disks: "local-disk ${disk_size} HDD"
  }
}

task check_is_binary {

  input {
    File pheno_file
    File pheno_list
  }

  Int disk_size = ceil(size(pheno_file,'GB')) + 2 * 2
  String outfile = "is_binary.tsv"

  command <<<
    set -eu
    PHENOLIST=$(tr -s ' \t' '\n' < ~{pheno_list} | awk 'NF && !seen[$0]++' | paste -sd,)

    # no pipefail: awk exits early once every pheno is resolved, which SIGPIPEs the still-writing zcat -- harmless, and
    # the sanity check below catches any case where a pheno didn't actually get classified.
    zcat -f ~{pheno_file} | awk -F'\t' -v phenolist="$PHENOLIST" '
      BEGIN { n=split(phenolist,names,","); remaining=n }
      NR==1 {
        for (i=1;i<=NF;i++) idx[$i]=i
        for (j=1;j<=n;j++) if (!(names[j] in idx)) { missing[names[j]]=1; remaining-- }
        if (remaining==0) exit
        next
      }
      {
        for (j=1;j<=n;j++) {
          p=names[j]
          if (p in missing || p in bad) continue
          v=$(idx[p])
          if (v=="" || v=="NA") continue
          if (v!=0 && v!=1) { bad[p]=1; remaining-- }
        }
        if (remaining==0) exit
      }
      END {
        for (j=1;j<=n;j++) {
          p=names[j]
          print p"\t"(p in missing ? "NA" : (p in bad ? 0 : 1))
        }
      }
    ' > ~{outfile}

    # every requested pheno must resolve to a real 0/1 classification
    while IFS=$'\t' read -r p v; do
      if [[ "$v" != "0" && "$v" != "1" ]]; then
        echo "ERROR: could not classify phenotype '$p' as binary/quantitative (got '$v')" >&2
        exit 1
      fi
    done < ~{outfile}

    cut -f2 ~{outfile} | sort | uniq -c | sort -nr >&2
  >>>

  output {
    File is_binary_tsv = outfile
  }

  runtime {
    disks: "local-disk ${disk_size} HDD"
  }
}

task filter_covariates {

  input {
    File pheno_file
    Array[String] covariates
    File pheno_list
    Int threshold_cov_count
  }

  String outfile = "./pheno_cov_map_" + threshold_cov_count + ".txt"
  Int disk_size = ceil(size(pheno_file,'GB')) + 2 * 2

  command <<<
    set -euxo pipefail

    python3 <<CODE
    import pandas as pd
    import numpy as np

    #read in phenos as list of phenos regardless
    tot_phenos = []
    phenos_groups = []
    with open('~{pheno_list}') as i:
        for line in i:
            phenos = line.strip().split()
            phenos_groups.append(phenos)
            tot_phenos += phenos

    #read in phenos mapping all valid entries to 1 and NAs to 0
    pheno_df = pd.read_csv('~{pheno_file}',sep='\t',usecols=tot_phenos).notna().astype(int)
    print(pheno_df)
    # read in covariates getting absolute values (handles PCs)
    covariates = '~{sep="," covariates}'.split(',')
    cov_df = pd.read_csv('~{pheno_file}',sep='\t',usecols=covariates).abs()
    print(cov_df)

    # now for each pheno calculate product of each covariate with itself

    with open('~{outfile}','wt') as o,open('~{outfile}'.replace('.txt','.err.txt'),'wt') as tmp_err:
        for i,pheno_list in enumerate(phenos_groups):
            pheno_name = ','.join(pheno_list)
            # for each group of phenos (possibly a single one) multiply all covs and pheno column and count how many non 0 entries are there: this means that the entry has a valid pheno and a non null covariates.
            df = pd.DataFrame()
            for pheno in pheno_list:
                m = (cov_df.mul(pheno_df[pheno],0)>0).sum().to_frame(pheno)
                df = pd.concat([df,m],axis=1)

            print(f"{i+1}/{len(phenos_groups)} {pheno_name}")
            #If it's a group of phenos the min function will return the lowest count across all phenos
            tmp_df = df[pheno_list].min(axis=1)
            covs = tmp_df.index[tmp_df >= ~{threshold_cov_count}].tolist()
            missing_covs = [elem for elem in covariates if elem not in covs]
            if missing_covs: tmp_err.write(f"{pheno_name}\t{','.join(missing_covs)}\n")
            o.write(f"{pheno_name}\t{','.join(covs)}\n")
    CODE
  >>>

  output {
    File cov_pheno_map = outfile
    File cov_pheno_warning = sub(outfile,'.txt','.err.txt')
  }

  runtime {
    memory: "64 GB"
    disks: "local-disk ${disk_size} HDD"
  }
}

task extract_cond_regions {

  input {
    String pheno
    String sumstats_root
    String region_root
    String chr_col
    String pos_col
    String ref_col
    String alt_col
    String mlogp_col
    Array[String] chroms
    Boolean add_hla
    Float mlogp_threshold
  }

  File region = sub(region_root,"PHENO",pheno)
  File sumstats = sub(sumstats_root,"PHENO",pheno)
  File index = sumstats + ".tbi"
  Int disk_size = ceil(size(sumstats,'GB')) + ceil(size(region,'GB')) + 1

  String outfile = pheno + "_sig_hits.txt"

  command <<<
    set -euo pipefail

    PHENO="~{pheno}"
    SUMSTATS="~{sumstats}"
    PVAL_THRESHOLD="~{mlogp_threshold}"
    CHR_COL="~{chr_col}"
    POS_COL="~{pos_col}"
    REF_COL="~{ref_col}"
    ALT_COL="~{alt_col}"
    MLOGP_COL="~{mlogp_col}"
    CHROMS=(~{sep=" " chroms})
    OUT="."

    REGIONS="tmp.bed"
    cat ~{region} > "$REGIONS"
    ~{if add_hla then "echo -e '6\t29000000\t34000000' >> tmp.bed" else ""}
    echo "$(IFS=' '; echo "${CHROMS[*]}") chromosomes included."
    echo "> ${PHENO} <"

    REGION_TMP="${OUT}/${PHENO}_subset_regions.txt.region.bed"
    CHROM_LIST=$(IFS=,; echo "${CHROMS[*]}")

    #subset regions to include only selected chroms
    awk -v chrom_list="${CHROM_LIST}" '
      BEGIN { n = split(chrom_list, arr, ","); for (i=1;i<=n;i++) valid[arr[i]] = 1 }
      ($1 in valid)
    ' "$REGIONS" > "$REGION_TMP"

    TOT_REGIONS=$(wc -l < "$REGION_TMP")
    [[ "$TOT_REGIONS" -gt 0 ]] && echo "${TOT_REGIONS} regions to be filtered." || echo "WARNING: ${REGIONS} empty. No hits will be returned!" >&2

    echo -n "Processing data..."

    #get column indices
    HEADER=$(tabix -H "$SUMSTATS")
    read -r P M R A <<< "$(awk -F'\t' -v pos_col="$POS_COL" -v mlogp_col="$MLOGP_COL" -v ref_col="$REF_COL" -v alt_col="$ALT_COL" '
      { for(i=1;i<=NF;i++) idx[$i]=i; print idx[pos_col], idx[mlogp_col], idx[ref_col], idx[alt_col] }
    ' <<< "$HEADER")"

    OUT_FILE="${OUT}/${PHENO}_sig_hits.txt"
    rm -f "${OUT}/${PHENO}_sig_hits_"*.txt "$OUT_FILE"
    touch "$OUT_FILE"

    while read -r chrom start end; do
      # subset sumstats to each region, filter by threshold, sort by mlogp, and take the top hit (if any) -- one line per region
      top=$(tabix "$SUMSTATS" "${chrom}:${start}-${end}" \
        | awk -F'\t' -v m="$M" -v t="$PVAL_THRESHOLD" '$m>t' \
        | sort -t$'\t' -k"${M},${M}rn" | { head -1; cat >/dev/null; })
      [[ -z "$top" ]] && continue

      IFS=$'\t' read -r -a f <<< "$top"
      mlogp="${f[$((M-1))]}" pos="${f[$((P-1))]}" ref="${f[$((R-1))]}" alt="${f[$((A-1))]}"
      # fix chromosome label for variant ID
      label="$chrom"; [[ "$chrom" == "23" ]] && label="X"
      variant="chr${label}_${pos}_${ref}_${alt}"
      region="${chrom}:${start}-${end}"
      line="${PHENO}\t${chrom}\t${region}\t${variant}\t${mlogp}"
      echo -e "$line" >> "$OUT_FILE"
      echo -e "$line" >> "${OUT}/${PHENO}_sig_hits_${chrom}.txt"
    done < "$REGION_TMP"

    echo "done."
    echo "dumping results to ${OUT_FILE}"
  >>>

  output {
    Array[File] pheno_chrom_regions = glob("*sig_hits_*")
    File gw_sig_res = outfile
  }

  runtime {
    disks: "local-disk ${disk_size} HDD"
  }
}

task merge_regions {

  input {
    Array[File] hits
    Boolean test
  }

  String outfile = "regions.txt"

  command <<<
    # test mode: cap each pheno's discovered regions to 2 random rows, so with the upstream
    # 10-pheno test-mode cap (validate_regions) the whole run is capped at 20 regions total.
    while read f; do
      if ~{true="true" false="false" test}; then
        shuf "$f" | head -n 2 >> tmp.txt
      else
        cat "$f" >> tmp.txt
      fi
    done < ~{write_lines(hits)}
    cat tmp.txt > ~{outfile}
  >>>

  output {
    File regions = outfile
  }
}
