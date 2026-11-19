version development

workflow regenie_lof {

  input {
    File null_map
    Float max_maf
    Float info_filter
    Array[String] lof_list
    File chunk_manifest  # release bgen chunks: path, chrom, first pos, last pos (regenie-pipelines scripts/return_bgen_chunks_limits.sh)
    Int n_batches = 40   # max number of shards subsetting the chunks
    File cov_file        # phenotypes + covariates (FID IID ...)
    String bargs_binary  # regenie step 2 options for binary phenos (0/1/NA only)
    String bargs_quant   # regenie step 2 options for quantitative phenos
    String docker
    String prefix
    Int mlogp_filter
    Int gene_variants_min_count
  }

  # get lof variants
  call extract_variants {
    input:
    docker = docker,
    max_maf = max_maf,
    info_filter = info_filter,
    lof_list = lof_list,
    gene_variants_min_count = gene_variants_min_count
  }

  # bgen chunks holding lof variants, split into at most n_batches batches
  call select_chunks {
    input:
      docker = docker,
      lof_variants = extract_variants.lof_variants,
      chunk_manifest = chunk_manifest,
      n_batches = n_batches
  }

  # subset each batch of chunks to its lof variants
  # batches.tsv columns: batch id, comma-separated chunk bgens, comma-separated chunk .bgi indexes.
  # The comma -> newline + write_lines/read_lines round trip turns a column into an Array[File]
  # (same trick as regenie-pipelines conditional-analysis).
  scatter (batch in read_tsv(select_chunks.batches)) {
    Array[File] batch_bgens = read_lines(write_lines([sub(batch[1], ",", "\n")]))
    Array[File] batch_bgis = read_lines(write_lines([sub(batch[2], ",", "\n")]))
    call subset_chunks {
      input:
        docker = docker,
        batch = batch[0],
        bgens = batch_bgens,
        bgis = batch_bgis,
        chunk_variants = select_chunks.chunk_variants
    }
  }

  # concatenate the batches (in order) into the lof bgen
  call merge_bgen {
    input:
      docker = docker,
      bgens = subset_chunks.batch_bgen,
      sample = select_chunks.sample,
      chunk_variants = select_chunks.chunk_variants
  }

  # null_map: pheno <tab> step 1 LOCO predictions (binary and quantitative phenos)
  Array[Array[String]] pheno_data = read_tsv(null_map)

  # binary (only 0/1/NA) or quantitative, one pass over the cov file for all phenos
  call check_is_binary {
    input:
      docker = docker,
      pheno_file = cov_file,
      pheno_list = write_lines(transpose(pheno_data)[0])
  }
  Map[String,String] is_binary_map = read_map(check_is_binary.is_binary_tsv)

  scatter(p_data in pheno_data) {
    String pheno = p_data[0]
    Boolean is_binary = is_binary_map[pheno] == "1"
    call regenie {
      input :
        pheno = pheno,
        null = p_data[1],
        cov_file = cov_file,
        bargs = if is_binary then bargs_binary else bargs_quant,
        prefix = prefix,
        lof_variants = extract_variants.lof_variants,
        sets = extract_variants.sets,
        mask = extract_variants.mask,
        bins = max_maf,
        lof_bgen = merge_bgen.lof_bgen,
        mlogp_filter = mlogp_filter
    }
  }

  call merge_results {
    input:
      docker = docker,
      regenie_results = regenie.results,
      comp_files = regenie.sig_res,
      prefix = prefix,
      logs = regenie.log,
      sets = extract_variants.sets,
      max_maf = max_maf,
      info = info_filter,
      lof_list = lof_list,
      mlogp_filter = mlogp_filter
  }

}


task merge_results {

  input {
    String docker
    Array[File] regenie_results
    Array[File] comp_files
    Array[File] logs
    File sets
    # README STUFF
    String prefix
    File lof_template
    String max_maf
    Int mlogp_filter
    Array[String] lof_list
    String info
  }

  String res_file = prefix + "_lof.txt" # file with all hits
  String sig_file = prefix + "_sig.txt" # file with sig hits (and variant comparison)
  String log_file = prefix + "_lof.log" # file with merged logs
  String sql_file = prefix + "_lof.sql.txt" # file for sql import 
  String readme   = prefix + "_lof_readme" # output readme 
  String var_file = prefix + "_lof_variants.txt" # list of variants used
  String checksum = prefix + "_check.log" # check that all ran or not
  

  command <<<
  # filter results
  paste <(echo PHENO) <(zcat  ~{regenie_results[0]} | head -n2 | tail -n1 | tr ' ' '\t')   > ~{res_file} # write header
  while read f
  do pheno=$(basename $f .regenie.gz |sed 's/~{prefix}_lof_//g' )  && zcat  $f | sed -E 1,2d |  awk -v pheno="$pheno" '{print pheno" "$0}' |  tr ' ' '\t'   >> tmp.txt
  done 	< ~{write_lines(regenie_results)}
  
  sort -grk 13 tmp.txt | grep -vw TEST_FAIL  >> ~{res_file}
  cat tmp.txt | grep -w TEST_FAIL > ~{log_file}
  # merge logs
  while read f
  do paste <( grep -q  "Elapsed" $f && echo 1 || echo 0 ) <(grep "phenoColList"  $f |  awk '{print $2}') <(echo $f| sed 's/\/cromwell_root/gs:\//g')  >> ~{checksum} && cat $f >> ~{log_file}
  done < ~{write_lines(logs)}

  head -n1 ~{comp_files[0]} > ~{sig_file}
  cut -f1,4 ~{sets} > ~{var_file}
  
  while read f
  do

      cat $f | sed -E 1d | awk '{print $5-$6"\t"$0}'   >> sig_tmp.txt
  done <  ~{write_lines(comp_files)}
  sort -rgk 1 sig_tmp.txt | cut -f2- | sed -e 's/ /\t/g' >> ~{sig_file}

  touch ~{sql_file}
  #GZIP FILES
  bgzip ~{log_file} -c > ~{log_file}.gz
  bgzip ~{res_file} -c > ~{res_file}.gz
  
  python3 <<EOF
  import sys,os,gzip,re
  import numpy as np

  # read in gene dict
  with open('~{sets}') as i:
    gene_dict = {}
    for line in i:
        gene,*_,variants = line.strip().split()
        gene_dict[gene] = variants
        
  with open('~{sql_file}','wt') as o, open('~{res_file}') as i:
    next(i) # go past header
    rel = re.findall(r'\d+','~{prefix}')[0] #get release number from prefix
    for line in i:
        pheno,_,_,gene,_,_,_,n,_,beta,se,_,mlogp,_ = line.strip().split()
        try:
            mlogp=float(mlogp)
            if float(mlogp) > float('~{mlogp_filter}'):
                gene = gene.split(".Mask")[0]
                data = map(str,[rel,pheno,gene,gene_dict[gene],np.power(10,-float(mlogp)),beta,se,n])
                o.write('"' + '","'.join(data) + '"\n')  
  
        except:
            print(pheno,gene,mlogp)
            
  lof_list = f"[{'~{sep="," lof_list}'}]"
  tags = [("[PREFIX]",'~{prefix}'),("[N_GENES]",str(len(gene_dict))),("[MAF]",'~{max_maf}'),("[INFO]",'~{info}'),('[LOF_LIST]',lof_list),('[MLOGP]','~{mlogp_filter}')]
  with open('~{lof_template}') as i,open('~{readme}','wt') as o:
    for line in i:
        for tag in tags:
            line = line.replace(tag[0],tag[1])
        o.write(line)     
  EOF  

  >>>
  
  runtime {
    docker: "~{docker}"
    cpu: 1
    disks:  "local-disk 10 HDD"
    memory: "2 GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 1
    maxRetries: 2
  }
   
  output {
    File all_hits = res_file +'.gz'
    File sig_hits = sig_file
    File sql_hits = sql_file
    File variants = var_file
    File read     = readme 
    File log      = log_file + '.gz'
    File check    = checksum
  }  
}

task regenie{

  input {
    String docker
    File lof_bgen
    File cov_file
    String pheno
    File null
    File sets
    File mask
    File lof_variants
    String covariates
    String bargs
    String prefix
    Int cpus
    String bins
    String mask_type
    String sumstats_root
    Int mlogp_filter
  }

  Int disk_size = ceil(size(lof_bgen,"GB"))*2 + 10
  File lof_index =  lof_bgen + ".bgi"
  File lof_sample = lof_bgen + ".sample"
  Map [String,File] pheno_map = {pheno:null}

  Int mem = ceil(size(lof_bgen,"GB"))*4*cpus
  String regenie_results = prefix + "_lof_" + pheno + ".regenie.gz"

  File sumstats = sub(sumstats_root,"PHENO",pheno)
  File tabix = sumstats + ".tbi"
  String pheno_comp = prefix + "_" + pheno + "_lof_sig_res.txt"

  
  command <<<
  CPUS=$(grep -c ^processor /proc/cpuinfo)

  # REGENIE
  cat ~{write_map(pheno_map)} > pred.txt
  time regenie --step 2 --out ./~{prefix}_lof --threads $CPUS  ~{bargs}  --bgen ~{lof_bgen} --sample ~{lof_sample} --pred pred.txt    --phenoFile ~{cov_file} --covarFile ~{cov_file} --phenoColList ~{pheno} --covarColList ~{covariates} --aaf-bins ~{bins}  --build-mask ~{mask_type} --mask-def ~{mask} --set-list ~{sets} --anno-file ~{lof_variants}
  echo ~{cpus} ~{mem} $CPUS
  wc -l pred.txt

  ####### SIG HITS FILTER
  
  # KEEP ONLY RELEVANT HITS (MLOGP > FILTER)
  echo "HITS WITH MLOGP > THRESHOLD"
  TMP_HITS="tmp_hits.txt"
  paste <(zcat ~{regenie_results} | awk '$12>~{mlogp_filter}'  | cut -d " " -f 3 | cut -d . -f 1 | awk '{print $0"_GENESTRING"}') <(zcat ~{regenie_results} | awk '$12>~{mlogp_filter}' | cut -d " " -f 1,9,12) | sed -e 's/ /\t/g' > $TMP_HITS

  # GET GENE TO VAR MAPPING FOR ALL GENES WITH HTIS
  cut -f 1,2 ~{lof_variants} | awk '{print $2"_GENESTRING\t"$1}' | sort -k1 | grep -wf <(cut -f1 $TMP_HITS ) > sig_genes.txt
  wc -l sig_genes.txt $TMP_HITS
  
  # FOR SIG GENES GET DATA BY TABIXING FROM SUMSTATS
  echo "RELEVANT GENES WITH SIG MLOGP"
  TMP_VAR="tmp_sig_genes_variants.txt"
  rm -f $TMP_VAR && touch $TMP_VAR
  while read -a  line # for each gene that is a significant HIT tabix from sumstats data
  do
      GENE=${line[0]}
      CHROM=${line[2]}
      SUMSTATS_CHROM=${CHROM/X/23}  # sumstats use 23 for X
      POS=${line[3]}
      REF=${line[4]}
      ALT=${line[5]}
      paste <(echo -e $GENE"_GENESTRING\tchr"$CHROM"_"$POS"_"$REF"_"$ALT) <(tabix ~{sumstats} $SUMSTATS_CHROM:$POS-$POS | grep -w $REF | grep -w $ALT | cut -f 8,9,11) >> $TMP_VAR
  done < <(cat  sig_genes.txt |  sed 's/chr//g' | tr '_' '\t' )

  # THIS IS THE PART TO BE JOINED WITH THE REGENIE RESULT (I COLLAPSE ALL VARIANT DATA FOR EACH GENE INTO A SINGLE LINE)
  TMP_GENE_VAR_HIT="tmp_sig_hits_gene.txt"
  rm -f $TMP_GENE_VAR_HIT && touch $TMP_GENE_VAR_HIT
  while read GENE
  do
      TOP_MLOGP=$( grep -w $GENE $TMP_VAR|  cut -f3 | sort -gr  | head -n1 ) # grep gene and keep highest MLOGP
      GENE_DATA=$(cat $TMP_VAR | grep -w $GENE | cut -f2- | tr '\t' ',' | tr '\n' ','| sed 's/.$//' ) # combine all lines into one
      echo -e $GENE'\t'$TOP_MLOGP'\t'$GENE_DATA >> $TMP_GENE_VAR_HIT
  done < <(cut -f1 $TMP_VAR | sort | uniq)

  # JOIN REGENIE RESULT WITH GENE VARIANT DATA
  PHENO_GENE_HITS="results.txt"
  echo -e "PHENO\tGENE\tCHROM\tGENE_BETA\tGENE_MLOGP\tTOP_VAR_MLOGP\tVAR1,MLOGP1,BETA1,AF1,VAR2..." > $PHENO_GENE_HITS 
  PHENO=~{pheno}
  join -t $'\t' <(sort -k1 $TMP_HITS) <(sort -k1 $TMP_GENE_VAR_HIT) | sed 's/_GENESTRING//g' | awk '{print "PNAME\t"$0}' | sed -e "s/PNAME/${PHENO}/g"  >> $PHENO_GENE_HITS

  mv $PHENO_GENE_HITS ~{pheno_comp}

  >>>

    runtime {
    docker: "~{docker}"
    cpu: "~{cpus}"
    disks:  "local-disk ~{disk_size} HDD"
    memory: "~{mem} GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 1
    maxRetries: 2
  }
  
  output {
    File results = regenie_results
    File sig_res = pheno_comp
    File log = prefix + "_lof.log"
  }
}


task check_is_binary {
  # From regenie-pipelines wdl/conditional-analysis/regenie_conditional_analysis.wdl (command unchanged):
  # a pheno is binary (1) if its column only holds 0/1/NA, quantitative (0) otherwise.

  input {
    String docker
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
    docker: "~{docker}"
    cpu: "1"
    disks: "local-disk ~{disk_size} HDD"
    memory: "2 GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 2
    maxRetries: 2
  }
}

task select_chunks {
  # Assigns each lof variant to the release bgen chunk whose [first pos, last pos] contains it and splits
  # the needed chunks into at most n_batches contiguous batches. Batches follow manifest order
  # (chr1..chr22, chrX, by position), so concatenating them in order keeps the bgen sorted.

  input {
    String docker
    File lof_variants
    File chunk_manifest
    Int n_batches
  }

  command <<<
  lof_variants="~{lof_variants}"
  chunk_manifest="~{chunk_manifest}"
  n_batches=~{n_batches}

  set -euo pipefail
  export LC_ALL=C

  # 1. Variant -> chunk. The chunks call X "chrX" (IDs chrX_pos_ref_alt), as extract_variants does; chr23 IDs (older
  #    lof_variants files) are translated. Output, in manifest order: chunk bgen path <tab> bgen variant ID
  awk -F'\t' -v OFS='\t' '
    NR == FNR {
      c = $2; n[c]++; i = n[c]
      start[c,i] = $3; end[c,i] = $4; path[c,i] = $1; rank[c,i] = FNR
      next
    }
    {
      split($1, v, "_"); c = (v[1] == "chr23") ? "chrX" : v[1]; pos = v[2] + 0
      for (i = 1; i <= n[c]; i++) if (pos >= start[c,i] && pos <= end[c,i]) break
      if (i > n[c]) { print "ERROR: no chunk contains " $1 > "/dev/stderr"; failed = 1; next }
      print rank[c,i], pos, path[c,i], c "_" v[2] "_" v[3] "_" v[4]
    }
    END { if (failed) exit 1 }
  ' "$chunk_manifest" "$lof_variants" | sort -k1,1n -k2,2n | cut -f3,4 > chunk_variants.tsv
  cut -f1 chunk_variants.tsv | uniq > chunks.txt

  # 2. Split the needed chunks into min(n_batches, chunks) contiguous batches whose sizes differ by at most one.
  #    Output: batch id <tab> chunk bgens <tab> chunk .bgi indexes (comma-separated)
  awk -v OFS='\t' -v n_batches="$n_batches" -v n="$(wc -l < chunks.txt)" '
    BEGIN { nb = (n_batches < n) ? n_batches : n }
    {
      b = int((NR - 1) * nb / n)
      sep = (count[b]++ > 0) ? "," : ""
      bgens[b] = bgens[b] sep $0
      bgis[b] = bgis[b] sep $0 ".bgi"
      last = b
    }
    END { for (i = 0; i <= last; i++) print sprintf("%03d", i), bgens[i], bgis[i] }
  ' chunks.txt > batches.tsv

  # all chunks share the same samples: take the sample file of the first one
  echo "$(head -n1 chunks.txt).sample" > sample_path.txt

  echo "$(wc -l < chunk_variants.tsv) variants in $(wc -l < chunks.txt) of $(wc -l < "$chunk_manifest") chunks, $(wc -l < batches.tsv) batches" >&2
  >>>

  runtime {
    docker: "~{docker}"
    cpu: "1"
    disks: "local-disk 10 HDD"
    memory: "2 GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 2
    maxRetries: 2
  }

  output {
    File chunk_variants = "chunk_variants.tsv"
    File batches = "batches.tsv"
    String sample = read_string("sample_path.txt")
  }
}

task subset_chunks {
  # Subsets each chunk of one batch to its lof variants (bgenix, using the chunk .bgi index)
  # and concatenates them, in batch order, into one bgen.

  input {
    String docker
    String batch
    Array[File] bgens
    Array[File] bgis
    File chunk_variants
  }

  # estimate from the first chunk only: size() on every chunk/index means one storage call per file
  Int disk_size = ceil(size(bgens[0],"MB") * length(bgens) / 1024) + 10

  command <<<
  batch="~{batch}"
  bgens_file="~{write_lines(bgens)}"
  bgis_file="~{write_lines(bgis)}"
  chunk_variants="~{chunk_variants}"

  set -euo pipefail
  export LC_ALL=C

  # link each chunk next to its own index (bgenix looks for <bgen>.bgi)
  mkdir -p chunks subsets
  paste "$bgens_file" "$bgis_file" | while IFS=$'\t' read -r bgen bgi; do
    ln -s "$bgen" "chunks/$(basename "$bgen")"
    ln -s "$bgi" "chunks/$(basename "$bgen").bgi"
  done

  # subset each chunk to its variants (matched on the bgen rsid field)
  n_chunks=$(wc -l < "$bgens_file"); i=0
  while read -r bgen; do
    name=$(basename "$bgen"); i=$((i + 1))
    awk -F'\t' -v c="$name" '{ p = $1; sub(/.*\//, "", p) } p == c { print $2 }' "$chunk_variants" > ids.txt
    bgenix -g "chunks/$name" -incl-rsids ids.txt > "subsets/$name" 2> bgenix.log || { cat bgenix.log >&2; exit 1; }
    echo "batch $batch: $i/$n_chunks $name ($(wc -l < ids.txt) variants)" >&2
  done < "$bgens_file"

  cat-bgen -clobber -og "batch_${batch}.bgen" -g $(while read -r bgen; do echo "subsets/$(basename "$bgen")"; done < "$bgens_file")
  >>>

  runtime {
    docker: "~{docker}"
    cpu: "4"
    disks: "local-disk ~{disk_size} SSD"  # SSD: localizing the chunks is the slow part
    memory: "4 GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 2
    maxRetries: 2
  }

  output {
    File batch_bgen = "batch_~{batch}.bgen"
  }
}

task merge_bgen {
  # Concatenates the batch bgens (in order), indexes the result and checks it holds exactly the lof variants.

  input {
    String docker
    Array[File] bgens
    File sample
    File chunk_variants
  }

  # estimate from the first batch only (one size() call); x2: inputs + merged lof.bgen
  Int disk_size = ceil(size(bgens[0],"MB") * length(bgens) / 1024)*2 + 10

  command <<<
  bgens_file="~{write_lines(bgens)}"
  sample="~{sample}"
  chunk_variants="~{chunk_variants}"

  set -euo pipefail
  export LC_ALL=C

  cat-bgen -clobber -og lof.bgen -g $(cat "$bgens_file")
  bgenix -g lof.bgen -index -clobber
  cp "$sample" lof.bgen.sample

  # the merged bgen must contain exactly the requested variants
  bgenix -g lof.bgen -list 2> /dev/null | awk -F'\t' '!/^#/ && $1 != "alternate_ids" { print $2 }' | sort > found.txt
  cut -f2 "$chunk_variants" | sort > wanted.txt
  if ! cmp -s wanted.txt found.txt; then
    echo "ERROR: lof.bgen does not match the requested variants" >&2
    comm -23 wanted.txt found.txt | sed 's/^/  missing: /' >&2
    comm -13 wanted.txt found.txt | sed 's/^/  unexpected: /' >&2
    exit 1
  fi
  echo "lof.bgen: $(wc -l < found.txt) variants, all requested variants present" >&2
  >>>

  runtime {
    docker: "~{docker}"
    cpu: "1"
    disks: "local-disk ~{disk_size} HDD"
    memory: "4 GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 2
    maxRetries: 2
  }

  output {
    File lof_bgen = "lof.bgen"
    File lof_index = "lof.bgen.bgi"
    File lof_sample = "lof.bgen.sample"
  }
}

task extract_variants {
  # Command = scripts/extract_variants.sh: the parameter block maps the WDL inputs to the script's
  # variables, the rest is the script body. Keep the two in sync.

  input {
    String docker
    File annot_file
    Array[String] lof_list
    Float max_maf
    Float info_filter
    Int gene_variants_min_count
  }

  Int disk_size = ceil(size(annot_file,'GB'))*2 + 10

  command <<<
  annot_file="~{annot_file}"
  max_maf=~{max_maf}
  info_filter=~{info_filter}
  gene_variants_min_count=~{gene_variants_min_count}
  lof_list=$(paste -sd, "~{write_lines(lof_list)}")  # comma-separated

  set -euo pipefail
  export LC_ALL=C  # locale-independent sort order

  # 1. LoF variants passing the MAF and INFO filters.
  #    Columns are looked up by header name (fails if one is missing).
  #    MAF = min(AF, 1-AF), so variants where ALT is the major allele are handled too.
  #    Variants with missing AF or INFO are dropped.
  #    Output, sorted by variant: chrN_pos_ref_alt (chrX for 23) <tab> gene <tab> consequence
  zcat -f "$annot_file" | awk -F'\t' -v OFS='\t' \
      -v lof_list="$lof_list" -v max_maf="$max_maf" -v info_filter="$info_filter" '
    BEGIN { n = split(lof_list, l, ","); for (i = 1; i <= n; i++) is_lof[l[i]] = 1 }
    NR == 1 {
      for (i = 1; i <= NF; i++) col[$i] = i
      n = split("#variant AF INFO most_severe gene_most_severe", req, " ")
      for (i = 1; i <= n; i++) if (!(req[i] in col)) { print "missing column in annot_file: " req[i] > "/dev/stderr"; failed = 1; exit 1 }
      next
    }
    !($col["most_severe"] in is_lof) { next }
    {
      n_lof++
      af = $col["AF"]; info = $col["INFO"]
      if (af == "" || af == "NA" || info == "" || info == "NA") { n_na++; next }
      maf = (af + 0 > 0.5) ? 1 - af : af + 0
      if (maf <= 0 || maf >= max_maf + 0) { n_maf++; next }
      if (info + 0 <= info_filter + 0)     { n_info++; next }
      id = "chr" $col["#variant"]; gsub(":", "_", id); sub(/^chr23_/, "chrX_", id)  # X as in the bgen chunks: chrX_pos_ref_alt
      print id, $col["gene_most_severe"], $col["most_severe"]
    }
    END {
      if (!failed) printf "LoF variants: %d. Removed: missing AF/INFO %d, MAF %d, INFO %d. Kept: %d\n", n_lof, n_na, n_maf, n_info, n_lof - n_na - n_maf - n_info > "/dev/stderr"
    }
  ' | sort -k1,1 > filtered.tsv

  # 2. Keep only genes with >= gene_variants_min_count variants (first pass counts, second pass keeps).
  #    This is the regenie annotation file: variant <tab> gene <tab> consequence
  awk -F'\t' -v min_count="$gene_variants_min_count" '
    NR == FNR { n[$2]++; next }
    n[$2] >= min_count + 0
  ' filtered.tsv filtered.tsv > lof_variants.txt

  # 3. Regenie set list: gene <tab> chrom <tab> pos <tab> variant1,variant2,...
  #    chrom/pos are taken from the gene's first variant.
  awk -F'\t' -v OFS='\t' '
    !($2 in vars) { genes[++n] = $2; vars[$2] = $1; next }
                  { vars[$2] = vars[$2] "," $1 }
    END { for (i = 1; i <= n; i++) { g = genes[i]; split(vars[g], p, "_"); print g, p[1], p[2], vars[g] } }
  ' lof_variants.txt | sort -k1,1 > sets.tsv

  # 4. Regenie mask definition: a single mask with every consequence present.
  printf 'Mask1\t%s\n' "$(cut -f3 lof_variants.txt | sort -u | paste -sd, -)" > mask.txt

  echo "Genes: $(wc -l < sets.tsv), variants: $(wc -l < lof_variants.txt), mask: $(cat mask.txt)" >&2
  >>>
  runtime {
    docker: "~{docker}"
    cpu: "1"
    disks:   "local-disk ~{disk_size} HDD"
    memory: "2 GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 2
    maxRetries: 2
  }

  output {
    File lof_variants = "lof_variants.txt"
    File mask = "./mask.txt"
    File sets = "./sets.tsv"
  }
}
