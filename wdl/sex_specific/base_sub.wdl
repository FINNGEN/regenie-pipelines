task step2 {

    Array[String] phenos
    File cov_pheno
    String covariates
    String test
    # regenie accepts either no test specified (i.e. additive), or '--test dominant' or '--test recessive'. '--test additive' is an error.
    String test_cmd = if test == "additive" then "" else "--test "+ test
    Boolean is_binary
    File bgen
    File bgi = bgen + ".bgi"
    File sample = bgen + ".sample"
    String grm_bed
    String prefix = basename(grm_bed, ".bed")
    Array[File] loco
    Array[File] nulls
    Int bsize
    String options
    File females
    File males
    String dollar="$"
    String docker
    String cpuplatform

    command <<<
        set -ex
        cp ${males} males
        cp ${females} females

        #create null and predlists
        paste ${write_lines(phenos)} ${write_lines(nulls)} > firth_list
        paste ${write_lines(phenos)} ${write_lines(loco)} > loco_list
        n_cpu=`grep -c ^processor /proc/cpuinfo`


        # Check status of checkpoint with python script
        cat << "__EOF__" > script.py
        import sys
        import glob
        variant_file = sys.argv[1]
        analysis_type = sys.argv[2].lower()
        endpoints = sys.argv[3].split(",")
        prefix= sys.argv[4]


        NO_VARIANTS_MIN_VALUE=1000000000
        NO_VARIANTS_MAX_VALUE=-1000000000
        def read_processed_range(files):
            """
            Read in the processed range from files
            Returns an optional tuple of chromosome, first position, last position
            Returns a tuple of chromosome, first position, last position
            """
            chroms = []
            pos_mins = []
            pos_maxs = []
            for fname in files:
                chrom = ""
                pos_min = NO_VARIANTS_MIN_VALUE
                pos_max = NO_VARIANTS_MAX_VALUE
                with open(fname,"r") as f:
                    l = f.readline()
                    c = l.strip().split(" ")
                    # if header is malformed, then nothing is processed.
                    if not ( (len(c) == 18 and analysis_type == "true") or (len(c)==14 and analysis_type == "false")):
                        return None
                    for line in f:
                        c = line.strip().split(" ")
                        #if line is malformed, break
                        if not ( (len(c) == 18 and analysis_type == "true") or (len(c)==14 and analysis_type == "false")):
                            break
                        if chrom == "":
                            chrom = c[0]
                        pos = int(c[1])
                        pos_min = min(pos,pos_min)
                        pos_max = max(pos,pos_max)
                pos_mins.append(pos_min)
                pos_maxs.append(pos_max)
                chroms.append(chrom)
            pos_min = min(pos_mins)
            pos_max = min(pos_maxs)
            if not all([a == chroms[0] for a in chroms]):
                #odd, all files don't have same chromosome
                print("Not all files have same chromosome! This is very odd!",file=sys.stderr)
            if pos_min != NO_VARIANTS_MIN_VALUE and pos_max != NO_VARIANTS_MAX_VALUE:
                return (chroms[0],pos_min,pos_max)
            return None
            
        def get_remaining_range(processed_range, all_variants_file)->str:
            """
            Get range still to be processed for files
            Returns a tuple of chromosome, first position, last position
            """
            # read range in all variants file
            chrom = ""
            pos_min = NO_VARIANTS_MIN_VALUE
            pos_max = NO_VARIANTS_MAX_VALUE
            with open(all_variants_file,"r") as f:
                _ = f.readline()
                for line in f:
                    c = line.strip().split("\t")
                    if chrom == "":
                        chrom = c[2].replace("chr","")
                    pos = int(c[3])
                    pos_min = min(pos,pos_min)
                    pos_max = max(pos,pos_max)
            #the processed range always should start from the beginning of the all variants range
            #if processed range is chr1:a-b, and complete range is chr1:a-c, where b < c, then the remaining range is chr1:b-c
            #we will want to be inclusive (i.e. b as start instead of b+1), since it is possible that there might be e.g. a multiallelic variant
            # and only one or part of its variants were included.
            return (chrom,processed_range[2],pos_max) 

        complete_male_files = glob.glob("checkpoint_folder/*.males*.regenie.gz")
        incomplete_male_files = glob.glob("checkpoint_folder/*.males*.regenie")
        male_log = glob.glob(f"checkpoint_folder/{prefix}.sex_spec.males.log")
        complete_female_files = glob.glob("checkpoint_folder/*.females*.regenie.gz")
        incomplete_female_files = glob.glob("checkpoint_folder/*.females*.regenie")
        female_log = glob.glob(f"checkpoint_folder/{prefix}.sex_spec.females")
        #check status of male and female endpoints
        males_run = all([any([f"{a}.regenie" in b for b in complete_male_files]) for a in endpoints]) and (len(male_log)==1)
        females_run = all([any([f"{a}.regenie" in b for b in complete_female_files]) for a in endpoints]) and (len(female_log)==1)
        if males_run:
            with open("males_done","w") as f:
                f.write("True")
        elif len(incomplete_male_files) != 0:
            # get processed range
            proc_range = read_processed_range(incomplete_male_files)
            if proc_range != None:
                remaining_range = get_remaining_range(proc_range,variant_file)
                with open("males_remaining_range","w") as f:
                    f.write(f"{remaining_range[0]}:{remaining_range[1]}-{remaining_range[2]}")
                with open("males_processed_range","w") as f:
                    f.write(f"{proc_range[0]}:{proc_range[1]}-{proc_range[2]}")
        if females_run:
            with open("females_done","w") as f:
                f.write("True")
        elif len(incomplete_female_files) != 0:
            proc_range = read_processed_range(incomplete_female_files)
            if proc_range != None:
                remaining_range = get_remaining_range(proc_range,variant_file)
                with open("females_remaining_range","a") as f:
                    f.write(f"{remaining_range[0]}:{remaining_range[1]}-{remaining_range[2]}")
                with open("females_processed_range","w") as f:
                    f.write(f"{proc_range[0]}:{proc_range[1]}-{proc_range[2]}")
        __EOF__

        # combine outputs from possibly multiple preempted runs
        cat << "__EOF__" > combine_outputs.py
        import sys
        import glob

        glob_command = sys.argv[1]
        output_name = sys.argv[2]
        phenos = sys.argv[3].split(",")
        files = glob.glob(glob_command)
        #now, read all files into memory
        data = []
        headers = []
        #do endpoints separately
        for p in phenos:
            files_ = [a for a in files if f"{p}.regenie" in a]
            for fname in files_:
                with open(fname,"r") as f:
                    headers.append(f.readline().strip().split(" "))
                    for l in f:
                        data.append(l.strip().split(" "))
            # check headers are the same
            if not all([a == headers[0] for a in headers]):
                print("NOT ALL HEADERS SAME!")
                for i in len(headers):
                    print(files_[i],header[i])
                sys.exit(1)
            header = headers[0]
            # order data by chrom:pos:ref:alt
            data2 = sorted(data,key=lambda x:(x[0],int(x[1])))
            #make unique
            data = []
            dataset = set()
            var = lambda x:(x[0],x[1],x[3],x[4])
            for d in data2:
                if var(d) in dataset:
                    continue
                else:
                    dataset.add(var(d))
                    data.append(d)
            with open(f"{output_name}.{p}","w") as of:
                of.write(" ".join(header)+"\n")
                for l in data:
                    of.write(" ".join(l)+"\n")
        __EOF__


        # write checkpointing script
        cat << "__EOF__" > checkpoint.sh
        #!/bin/bash
        touch placeholder.regenie.placeholder placeholder.log
        while [ 1 -eq 1 ];do
            tar -cf checkpoint_new.tar *.regenie* *.log
            cp checkpoint.tar checkpoint_old.tar
            mv checkpoint_new.tar checkpoint.tar
            sleep 60
        done
        __EOF__
        chmod +x checkpoint.sh
        # TODO: get range of variants
        # get list of variants with bgenix
        bgenix -g ${bgen} -list |grep -Ev "^#" > variants.list
        # TODO: check status of checkpoint
        #create "workspace" for checkpoint
        mkdir checkpoint_folder
        # if checkpoint is available, untar it
        if test -f checkpoint.tar;then
            if tar -xf checkpoint.tar -C checkpoint_folder/; then
                python3 script.py variants.list "${is_binary}" "${sep="," phenos}" "${prefix}"
            fi
        fi

        # start checkpointing script
        ./checkpoint.sh &
        CHECKPOINT_PID=$!

        # check if checkpoint file contains male regenie outputs
        if test -f males_done; then
            #copy files to here
            cp checkpoint_folder/*.males*.regenie.gz ./
            cp checkpoint_folder/${prefix}.sex_spec.males.log ./
        else
            # test if we have processed some range
            if test -f males_remaining_range; then
                #move all non-complete male regenie outputs to have the range that they
                PROCESSED_M_RANGE=$(cat males_processed_range)
                for f in checkpoint_folder/*.males*.regenie;do
                    cp $f $f.$PROCESSED_M_RANGE
                done
                #copy all files with any range from checkpoint folder to main folder
                # this will fix the problem with one or more ranges missing when  
                cp checkpoint_folder/*.males_*regenie.* ./
                cat checkpoint_folder/*.males.log* > ./${prefix}.sex_spec.males.log.$PROCESSED_M_RANGE
                REMAINING_M_RANGE=$(cat males_remaining_range)
                #we have, use the calculated remaining range 
                M_RANGE_OPTION="--range "$REMAINING_M_RANGE
                echo "Remaining male range "$REMAINING_M_RANGE
            fi
            echo "running males analysis"
            echo "$(wc -l males) individuals in file"
            regenie \
            --step 2 \
            ${test_cmd} \
            ${if is_binary then "--bt --af-cc" else ""} \
            --bgen ${bgen} \
            --ref-first \
            --sample ${sample} \
            --keep males \
            --covarFile ${cov_pheno} \
            --covarColList ${covariates} \
            --phenoFile ${cov_pheno} \
            --phenoColList ${sep="," phenos} \
            --pred loco_list \
            ${if is_binary then "--use-null-firth firth_list" else ""} \
            --bsize ${bsize} \
            --threads $n_cpu \
            --out ${prefix}".sex_spec.males" \
            $M_RANGE_OPTION \
            ${options}
            #update checkpoint: There are no males done, and therefore no females either.
            #no compression, as the biggest source of data is already compressed
            # stitch files together if necessary
            python3 combine_outputs.py "${prefix}.sex_spec.males_*.regenie*" temp_males.regenie "${sep="," phenos}"
            for p in ${sep=" " phenos};do
                cat temp_males.regenie.$p |gzip >  ${prefix}.sex_spec.males_$p.regenie.gz
            done
            #compress
            tar -cf checkpoint.tar *.males*.regenie.gz ${prefix}.sex_spec.males.log*
        fi 

        # check if checkpoint contains female regenie outputs
        if test -f females_done; then
            #copy files to here. Males have already been copied, or they already exist in this folder.
            cp checkpoint_folder/*.females*.regenie.gz ./
            cp checkpoint_folder/${prefix}.sex_spec.females.log ./
        else
            # test if we have processed some range
            if test -f females_remaining_range; then
                #move all non-complete male regenie outputs to have the range that they
                PROCESSED_F_RANGE=$(cat females_processed_range)
                for f in checkpoint_folder/*.females*.regenie;do
                    cp $f $f.$PROCESSED_F_RANGE
                done
                cp checkpoint_folder/*.females_*regenie.* ./
                cat checkpoint_folder/*.females.log* > ./${prefix}.sex_spec.females.log.$PROCESSED_F_RANGE
                REMAINING_F_RANGE=$(cat females_remaining_range)
                #we have, use the calculated remaining range 
                F_RANGE_OPTION="--range "$REMAINING_F_RANGE
                echo "Remaining female range "$REMAINING_M_RANGE
            fi
            echo "running females analysis"
            echo "$(wc -l females) individuals in file"
            regenie \
            --step 2 \
            ${test_cmd} \
            ${if is_binary then "--bt --af-cc" else ""} \
            --bgen ${bgen} \
            --ref-first \
            --sample ${sample} \
            --keep females \
            --covarFile ${cov_pheno} \
            --covarColList ${covariates} \
            --phenoFile ${cov_pheno} \
            --phenoColList ${sep="," phenos} \
            --pred loco_list \
            ${if is_binary then "--use-null-firth firth_list" else ""} \
            --bsize ${bsize} \
            --threads $n_cpu \
            --out ${prefix}".sex_spec.females" \
            $F_RANGE_OPTION \
            ${options}
            #update checkpoint: All files have been done.
            python3 combine_outputs.py "${prefix}.sex_spec.females_*.regenie*" temp_females.regenie "${sep="," phenos}"
            for p in ${sep=" " phenos};do
                cat temp_females.regenie.$p |gzip >  ${prefix}.sex_spec.females_$p.regenie.gz
            done
            #no compression, as the biggest source of data is already compressed
            tar -cf checkpoint.tar *.males*.regenie.gz *.females*.regenie.gz ${prefix}.sex_spec.males.log* ${prefix}.sex_spec.females.log*
        fi 
        
        cat ${prefix}.sex_spec.males.log* ${prefix}.sex_spec.females.log* > ${prefix}.log 
        set +ex
        #stop checkpoint script
        kill $CHECKPOINT_PID

        echo "all files"
        ls * || true
        
    >>>

    output {
        File log = prefix + ".log"
        Array[File] male_files = glob("*.males*.regenie.gz")
        Array[File] female_files = glob("*.females*.regenie.gz")
    }

    runtime {
        docker: "${docker}"
        cpu: 4
        memory: "4 GB"
        disks: "local-disk " + (ceil(size(bgen, "G")) + 5) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 4
        noAddress: true
        #cpuPlatform:cpuplatform
        checkpointFile: "checkpoint.tar"

    }
}

task gather {

    Boolean is_binary
    Array[File] files_male
    Array[File] files_female

    String docker
    String cpuplatform
    command <<<

        set -euxo pipefail

        pheno=`basename ${files_male[0]} .regenie.gz | awk -F "." '{sub(/[^_]*_/, "", $NF); print $NF}'`
        mkdir regenie munged

        echo -e "`date`\tconcatenating result pieces into regenie/$pheno.regenie.gz, sorting by chr pos just in case"
        cat \
        <(zcat ${files_male[0]} | head -1) \
        <(for file in ${sep=" " files_male}; do
            zcat $file | tail -n+2
        done | sort -k1,1g -k2,2g) | bgzip > regenie/$pheno.males.regenie.gz

        cat \
        <(zcat ${files_female[0]} | head -1) \
        <(for file in ${sep=" " files_female}; do
            zcat $file | tail -n+2
        done | sort -k1,1g -k2,2g) | bgzip > regenie/$pheno.females.regenie.gz

        if [[ "${is_binary}" == "true" ]]; then
            echo -e "`date`\tconverting to munged/$pheno.males.gz to a format used for importing to pheweb, omitting variants with -log10p NA (unsuccessful Firth/SPA)"
            zcat regenie/$pheno.males.regenie.gz | awk '
            BEGIN {FS=" "; OFS="\t"; split("CHROM GENPOS ALLELE0 ALLELE1 LOG10P BETA SE A1FREQ A1FREQ_CASES A1FREQ_CONTROLS", REQUIRED_FIELDS)}
            NR==1 {for(i=1;i<=NF;i++) h[$i]=i;
                   for(i in REQUIRED_FIELDS) if (!(REQUIRED_FIELDS[i] in h)) {print REQUIRED_FIELDS[i]" expected in regenie header">>"/dev/stderr"; exit 1}
                   print "#chrom","pos","ref","alt","pval","mlogp","beta","sebeta","af_alt","af_alt_cases","af_alt_controls"}
            NR>1 && $h["LOG10P"]!="NA" {print $h["CHROM"],$h["GENPOS"],$h["ALLELE0"],$h["ALLELE1"],10^-$h["LOG10P"],$h["LOG10P"],$h["BETA"],$h["SE"],$h["A1FREQ"],$h["A1FREQ_CASES"],$h["A1FREQ_CONTROLS"]}' \
            | bgzip > munged/$pheno.males.gz

            zcat regenie/$pheno.females.regenie.gz | awk '
            BEGIN {FS=" "; OFS="\t"; split("CHROM GENPOS ALLELE0 ALLELE1 LOG10P BETA SE A1FREQ A1FREQ_CASES A1FREQ_CONTROLS", REQUIRED_FIELDS)}
            NR==1 {for(i=1;i<=NF;i++) h[$i]=i;
                   for(i in REQUIRED_FIELDS) if (!(REQUIRED_FIELDS[i] in h)) {print REQUIRED_FIELDS[i]" expected in regenie header">>"/dev/stderr"; exit 1}
                   print "#chrom","pos","ref","alt","pval","mlogp","beta","sebeta","af_alt","af_alt_cases","af_alt_controls"}
            NR>1 && $h["LOG10P"]!="NA" {print $h["CHROM"],$h["GENPOS"],$h["ALLELE0"],$h["ALLELE1"],10^-$h["LOG10P"],$h["LOG10P"],$h["BETA"],$h["SE"],$h["A1FREQ"],$h["A1FREQ_CASES"],$h["A1FREQ_CONTROLS"]}' \
            | bgzip > munged/$pheno.females.gz
        else
            echo -e "`date`\tconverting to munged/$pheno.males.gz to a format used for importing to pheweb"
            zcat regenie/$pheno.males.regenie.gz | awk '
            BEGIN {FS=" "; OFS="\t"; split("CHROM GENPOS ALLELE0 ALLELE1 LOG10P BETA SE A1FREQ", REQUIRED_FIELDS)}
            NR==1 {for(i=1;i<=NF;i++) h[$i]=i;
                   for(i in REQUIRED_FIELDS) if (!(REQUIRED_FIELDS[i] in h)) {print REQUIRED_FIELDS[i]" expected in regenie header">>"/dev/stderr"; exit 1}
                   print "#chrom","pos","ref","alt","pval","mlogp","beta","sebeta","af_alt"}
            NR>1  {print $h["CHROM"],$h["GENPOS"],$h["ALLELE0"],$h["ALLELE1"],10^-$h["LOG10P"],$h["LOG10P"],$h["BETA"],$h["SE"],$h["A1FREQ"]}' \
            | bgzip > munged/$pheno.males.gz

            echo -e "`date`\tconverting to munged/$pheno.males.gz to a format used for importing to pheweb"
            zcat regenie/$pheno.females.regenie.gz | awk '
            BEGIN {FS=" "; OFS="\t"; split("CHROM GENPOS ALLELE0 ALLELE1 LOG10P BETA SE A1FREQ", REQUIRED_FIELDS)}
            NR==1 {for(i=1;i<=NF;i++) h[$i]=i;
                   for(i in REQUIRED_FIELDS) if (!(REQUIRED_FIELDS[i] in h)) {print REQUIRED_FIELDS[i]" expected in regenie header">>"/dev/stderr"; exit 1}
                   print "#chrom","pos","ref","alt","pval","mlogp","beta","sebeta","af_alt"}
            NR>1  {print $h["CHROM"],$h["GENPOS"],$h["ALLELE0"],$h["ALLELE1"],10^-$h["LOG10P"],$h["LOG10P"],$h["BETA"],$h["SE"],$h["A1FREQ"]}' \
            | bgzip > munged/$pheno.females.gz
        fi

        echo -e "`date`\tprinting only chr, pos, pval to speed up and reduce memory use of qqplot.R"
        zcat munged/$pheno.males.gz | cut -f1,2,5 > $pheno.males
        zcat munged/$pheno.females.gz | cut -f1,2,5 > $pheno.females

        echo -e "`date`\trunning qqplot.R"
        qqplot.R --file $pheno.males --chrcol "#chrom" --bp_col "pos" --pval_col "pval" --loglog_pval 10 > qq_out_males 2> qq_err_males
        qqplot.R --file $pheno.females --chrcol "#chrom" --bp_col "pos" --pval_col "pval" --loglog_pval 10 > qq_out_females 2> qq_err_females

        echo -e "`date`\ttabixing munged/$pheno.males.gz"
        tabix -s1 -b2 -e2 munged/$pheno.males.gz
        tabix -s1 -b2 -e2 munged/$pheno.females.gz
        echo -e "`date`\tdone"

    >>>

    output {
        Array[File] regenie = glob("regenie/*.regenie.gz")
        Array[File] pheweb = glob("munged/*.gz")
        Array[File] pheweb_tbi = glob("munged/*.gz.tbi")
        Array[File] qq_out = glob("qq_out*")
        Array[File] qq_err = glob("qq_err*")
        Array[File] pngs = glob("*.png")
        Array[File] quantiles = glob("*qquantiles.txt")
    }

    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "6 GB"
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
        #cpuPlatform:cpuplatform
    }
}

task summary{

    Array[File] input_files
    File finngen_annotation
    File gnomad_annotation
    File finngen_tbi = finngen_annotation + ".tbi"
    File gnomad_tbi = gnomad_annotation + ".tbi"

    Float summary_pval_thresh
    Float coding_pval_thresh
    String docker

    command <<<
        set -euxo pipefail

        cat << "__EOF__" > script.py
        import sys
        #read file
        fname = sys.argv[1]
        finngen_annotation_file = sys.argv[2]
        gnomad_annotation_file  = sys.argv[3]
        summary_threshold = ${summary_pval_thresh}
        coding_threshold = ${coding_pval_thresh}

        import pysam
        import gzip
        import os
        from typing import NamedTuple

        coding_groups = set([
            "transcript_ablation",
            "splice_donor_variant",
            "stop_gained",
            "splice_acceptor_variant",
            "frameshift_variant",
            "stop_lost",
            "start_lost",
            "inframe_insertion",
            "inframe_deletion",
            "missense_variant",
            "protein_altering_variant"
        ])

        pheno = os.path.splitext(os.path.basename(fname))[0]
        output_summary_name = pheno + "_summary.txt"
        output_coding_name = pheno + "_coding.txt"
        FGAnnotation = NamedTuple('FGAnnotation',[('gene',str),('consequence',str)])
        GnomadAnnotation = NamedTuple('GnomadAnnotation',[('finaf',str),('nfeaf',str),('enrich',str),('rsid',str)])

        def get_header(reader,file):
            with reader(file, "rt") as f:
                l = f.readline()
                return l.strip("\n").split("\t")

        def get_fg_annotation(iterator,variant,gene_idx, consequence_idx,variant_idx) -> FGAnnotation:
            for v in iterator:
                cols = v.strip("\n").split("\t")
                if cols[variant_idx] == variant:
                    return FGAnnotation(cols[gene_idx],cols[consequence_idx])
            return FGAnnotation("","")

        def get_gnomad_annotation(iterator,cpra,finaf_idx, nfeaf_idx, rsid_idx,gd_cpra) -> GnomadAnnotation:
            for v in iterator:
                cols = v.strip("\n").split("\t")

                if (cols[gd_cpra[0]] == "chr"+cpra[0]) and (cols[gd_cpra[1]] == str(cpra[1])) and (cols[gd_cpra[2]] == cpra[2]) and (cols[gd_cpra[3]] == cpra[3]):
                    return GnomadAnnotation(cols[finaf_idx],cols[nfeaf_idx],cols[enrich_idx],cols[rsid_idx])
            return GnomadAnnotation(".",".",".",".")

        #required columns
        fg_req_cols=["#variant","gene_most_severe","most_severe"]
        gd_req_cols=["fin.AF","nfsee.AF","enrichment_nfsee","rsid"]
        #open finngen annotation tabix
        fg_tabix = pysam.TabixFile(finngen_annotation_file,parser=None)
        #get fg header column positions
        fg_header = get_header(gzip.open, finngen_annotation_file)
        fg_idx = {v:i for (i,v) in enumerate(fg_header)}
        #open gnomad annotation tabix
        gnomad_tabix = pysam.TabixFile(gnomad_annotation_file,parser=None)
        gd_header = get_header(gzip.open, gnomad_annotation_file)
        gd_idx = {v:i for (i,v) in enumerate(gd_header)}
        gd_cpra = (
            gd_idx["chrom"],
            gd_idx["pos"],
            gd_idx["ref"],
            gd_idx["alt"]
        )
        #check for fg column existence
        if not all([a in fg_header for a in fg_req_cols]):
            raise Exception("Not all columns present in FinnGen annotation! Aborting...")
        #check for gnomad column existence
        if not all([a in gd_header for a in gd_req_cols]):
            raise Exception("Not all columns present in Gnomad annotation! Aborting...")
        var_idx, gene_idx, cons_idx = (fg_idx["#variant"],fg_idx["gene_most_severe"],fg_idx["most_severe"])

        finaf_idx, nfeaf_idx, enrich_idx, rsid_idx = (gd_idx["fin.AF"],gd_idx["nfsee.AF"],gd_idx["enrichment_nfsee"],gd_idx["rsid"])



        with gzip.open(fname, "rt") as file:
            #open output file
            with open(output_summary_name,"w") as summary_outfile, open(output_coding_name,"w") as coding_outfile:
                #read header
                header = file.readline().strip("\n").split('\t')
                #find p-value index
                header_idx = {v:i for (i,v) in enumerate(header)}
                pval_idx = header_idx["pval"]
                cid,pid,rid,aid = (
                    header_idx["#chrom"],
                    header_idx["pos"],
                    header_idx["ref"],
                    header_idx["alt"]
                )

                #add gene name, consequence, gnomad finnish af, nfsee af, rsid
                header.extend(["fin.AF","nfsee.AF","rsid","gene_most_severe","most_severe"])
                summary_outfile.write("\t".join(header)+"\n")
                #add pheno, fin enrichment
                header.extend(["fin.enrichment","phenotype"])
                coding_outfile.write("\t".join(header)+"\n")
                
                #read lines
                for line in file:
                    line_columns = line.strip("\n").split('\t')
                    pvalue = float(line_columns[pval_idx])
                    if pvalue < coding_threshold or pvalue < summary_threshold:
                        cpra= (line_columns[cid],int(float(line_columns[pid])),line_columns[rid],line_columns[aid])
                        variant = "{}:{}:{}:{}".format(cpra[0],cpra[1],cpra[2],cpra[3])
                        fg_c = cpra[0].replace("chr","").replace("X","23").replace("Y","24").replace("M","25").replace("MT","25")
                        gd_c = "chr" + cpra[0].replace("chr","").replace("23","X").replace("24","Y").replace("25","M")
                        #annotate
                        fg_iter = fg_tabix.fetch(fg_c,cpra[1]-1, cpra[1])
                        fg_a = get_fg_annotation(fg_iter,variant, gene_idx,cons_idx,var_idx)

                        #annotate
                        gnomad_iter = gnomad_tabix.fetch(gd_c,cpra[1]-1, cpra[1])
                        gd_a = get_gnomad_annotation(gnomad_iter,cpra,finaf_idx, nfeaf_idx, rsid_idx,gd_cpra)

                        line_columns.extend([
                            gd_a.finaf,
                            gd_a.nfeaf,
                            gd_a.rsid,
                            fg_a.gene,
                            fg_a.consequence,
                        ])
                        #gather row
                        #write to file
                        if pvalue < summary_threshold:
                            summary_outfile.write("\t".join(line_columns)+"\n")

                        #coding out
                        if pvalue < coding_threshold and fg_a.consequence in coding_groups:
                            line_columns.extend([
                                gd_a.enrich,
                                pheno
                            ])
                            coding_outfile.write("\t".join(line_columns)+"\n")
        print("summary created")
        __EOF__

        for phenofile in ${sep=" " input_files};do
            python3 script.py $phenofile ${finngen_annotation} ${gnomad_annotation}
        done

    >>>

    output {
        Array[File] summary_out = glob("*_summary.txt")
        Array[File] coding_out = glob("*_coding.txt")
    }

    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
    }
}
workflow regenie_step2 {

    String docker
    Array[String] phenos
    Boolean is_binary
    String cov_pheno
    String covariates

    Array[String] loco
    Array[String] nulls
    File females
    File males

    File bgenlist
    Array[String] bgens = read_lines(bgenlist)
    String cpuplatform
    scatter (bgen in bgens) {
        call step2 {
            input: docker=docker, phenos=phenos, is_binary=is_binary, cov_pheno=cov_pheno, covariates=covariates, loco=loco, nulls=nulls,bgen=bgen,females=females,males=males,cpuplatform=cpuplatform
        }
    }

    Array[Array[String]] pheno_results_males = transpose(step2.male_files)
    Array[Array[String]] pheno_results_females = transpose(step2.female_files)
    scatter (i in range(length(pheno_results_males))) {
        call gather {
            input: is_binary=is_binary, files_male=pheno_results_males[i],files_female=pheno_results_females[i],cpuplatform=cpuplatform
        }
        call summary {
            input: input_files=gather.pheweb
        }
    }

    output {
        Array[File] pheweb = flatten(gather.pheweb)
        Array[File] pheweb_tbi = flatten(gather.pheweb_tbi)
        Array[File] qq_out = flatten(gather.qq_out)
        Array[File] qq_err = flatten(gather.qq_err)
        Array[File] pngs = flatten(gather.pngs)
        Array[Array[File]] quantiles = gather.quantiles
        Array[File] summaries = flatten(summary.summary_out)
        Array[File] coding = flatten(summary.coding_out)
    }
}
