task step2 {

    Array[String] phenolist
    File cov_pheno
    String covariates
    String test
    # regenie accepts either no test specified (i.e. additive), or '--test dominant' or '--test recessive'. '--test additive' is an error.
    String test_cmd = if test == "additive" then "" else "--test "+ test
    Boolean is_binary
    File bgen
    File bgi = bgen + ".bgi"
    File sample = bgen + ".sample"
    File pred
    File firth_list
    String prefix = sub(basename(pred), "_pred.list", "") + "." + basename(bgen)
    Array[File] loco
    Array[File] nulls
    Int bsize
    String options

    String docker


    Float sex_specific_logpval
    String sex_col_name
    Boolean run_sex_specific

    String DOLLAR="$"

    command <<<

    ## create inputs
    n_cpu=`grep -c ^processor /proc/cpuinfo`
    # create loco and firth lists
    paste ${write_lines(phenolist)} ${write_lines(nulls)} > firth_list
    paste ${write_lines(phenolist)} ${write_lines(loco)} > loco_list
    ## write scripts into files

# Check status of checkpoint with python script
cat << "__EOF__" > checkpoint_status.py
import sys
import glob
import argparse
from enum import Flag,auto
from typing import Optional, Tuple, List
# chromosomal range type alias for convenience
chromRange = Tuple[str,int,int]
class AnalysisState(Flag):
    NOT_STARTED = auto()
    INCOMPLETE = auto()
    FINISHED = auto()

class SexSpecState(Flag):
    NOT_COMPLETE = auto()
    FINISHED = auto()


parser = argparse.ArgumentParser()
parser.add_argument("variant_list")
parser.add_argument("--analysis-type")
parser.add_argument("--endpoints")
parser.add_argument("--prefix")
args = parser.parse_args()
variant_file = args.variant_list
analysis_type = args.analysis_type.lower()
endpoints = args.endpoints.split(",")
prefix= args.prefix


NO_VARIANTS_MIN_VALUE=1000000000
NO_VARIANTS_MAX_VALUE=-1000000000
def read_processed_range(files:List[str])->Optional[chromRange]:
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
        print(f"Processed chroms:{chroms}")
    if pos_min != NO_VARIANTS_MIN_VALUE and pos_max != NO_VARIANTS_MAX_VALUE:
        return (chroms[0],pos_min,pos_max)
    return None
    
def get_chunk_range(all_variants_file:str)->chromRange:
    """
    Get range still to be processed for files
    Returns a tuple of chromosome, first position, last position
    """
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
    return (chrom,pos_min,pos_max) 

complete_files = glob.glob("checkpoint_folder/*.regenie.gz")
incomplete_files = glob.glob("checkpoint_folder/*.regenie")
# sex specific
sex_specific_files = glob.glob("checkpoint_folder/*.sex_spec.gz")
log = glob.glob(f"checkpoint_folder/{prefix}.log")
# all endpoints are run, if all endpoints are in complete files (gzipped files), and if there is a log.
endpoints_run :bool= all([any([f"{a}.regenie" in b for b in complete_files]) for a in endpoints]) and (len(log)==1)
incomplete_files_for_each_endpoint :bool = all([any([f"{a}.regenie" in b for b in incomplete_files]) for a in endpoints])
# determine current state
analysis_state:AnalysisState
if endpoints_run:
    analysis_state = AnalysisState.FINISHED
elif incomplete_files_for_each_endpoint:
    analysis_state = AnalysisState.INCOMPLETE
else:
    analysis_state = AnalysisState.NOT_STARTED

print(f"The state of analysis in checkpoint is {analysis_state}",file=sys.stderr)

if analysis_state == AnalysisState.FINISHED:
    with open("phenos_done","w") as f:
        f.write("True\n")
# if there are incomplete files for each of the endpoints
elif analysis_state == AnalysisState.INCOMPLETE:
    # get processed range
    proc_range = read_processed_range(incomplete_files)
    if proc_range != None:
        chunk_range = get_chunk_range(variant_file)
        remaining_range = (chunk_range[0],proc_range[2],chunk_range[2])
        with open("remaining_range","w") as f:
            f.write(f"{remaining_range[0]}:{remaining_range[1]}-{remaining_range[2]}")
        with open("processed_range","w") as f:
            f.write(f"{proc_range[0]}:{proc_range[1]}-{proc_range[2]}")
else:
    pass
# sex-specific endpoints
# I only look at the actual finished files, no male or female files separately. Let that be the MVP.
# I can also just check whether the file is done in checkpointing folder, and if so, just copy it out from there.
# That would be done in the actual script.


done_sex_spec = [a for a in endpoints if any([f"{a}.sex_spec.gz" in b for b in sex_specific_files])]
not_done_sex_spec = [a for a in endpoints if a not in done_sex_spec]
all_sex_spec_done = all([a in done_sex_spec for a in endpoints])
sex_spec_state:SexSpecState
if all_sex_spec_done:
    sex_spec_state = SexSpecState.FINISHED
else:
    sex_spec_state = SexSpecState.NOT_COMPLETE
print(f"The state of sex-specific analysis in checkpoint is {sex_spec_state}",file=sys.stderr)
if sex_spec_state == SexSpecState.FINISHED:
    with open("sex_spec_done","w") as f:
        f.write("True\n")
elif sex_spec_state == SexSpecState.NOT_COMPLETE:
    with open("not_done_sex_spec","w",encoding="utf-8") as f:
        for p in not_done_sex_spec:
            f.write(f"{p}\n")
__EOF__

# combine outputs from possibly multiple preempted runs
cat << "__EOF__" > combine_outputs.py
import sys
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("glob_command")
parser.add_argument("--output-name")
parser.add_argument("--endpoints")
args = parser.parse_args()
glob_command = args.glob_command
output_name = args.output_name
phenos = args.endpoints.split(",")

files = glob.glob(glob_command)
#do endpoints separately
for p in phenos:
    data = []
    headers = []
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
    # order data by chrom:pos
    data2 = sorted(data,key=lambda x:(x[0],int(x[1])))
    #make unique. Only include variants that are not already in variant set.
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

# write sex spec modifying script
echo '''import pandas as pd
import gzip
import sys
import math
from scipy.stats import norm
import numpy as np

pheno=sys.argv[1]
prefix=sys.argv[2]
base=sys.argv[3]

male=prefix+ ".sex_spec.males_" + pheno + ".gz"
female=prefix+".sex_spec.females_" + pheno + ".gz"


## NOTE:  See above where Echoing the expected header in case no sex specific results (above continue statement) and we want each shard to
## NOTE: have correct header so in gather step the first line of first file can serve as a header.
## NOTE: IF changing output columns in here, change the columns accordingly above.

sex_cols = ["ID","N","BETA","SE","LOG10P"]

basic = pd.read_csv( gzip.open(base), sep=" ")

cols = list(basic.columns)

binarycols = ["A1FREQ_CONTROLS","A1FREQ_CASES"] 

malestat = pd.read_csv( gzip.open( male ), sep=" " )
femalestat = pd.read_csv( gzip.open(female), sep=" ")

if  all ( [ c in cols for c in binarycols ]):
    ## add case control afs... 
    sex_cols.extend(binarycols)
    
malestat = malestat[sex_cols]
femalestat = femalestat[sex_cols]

combs = basic.merge(malestat.add_prefix("males_"),
          left_on="ID",right_on="males_ID").merge(femalestat.add_prefix("females_"),left_on="ID", right_on="females_ID")

combs["diff_beta"] = combs["males_BETA"]-combs["females_BETA"]

if(len(combs.index))>0:
  combs["p_diff"] = -((norm.logsf( abs(combs["diff_beta"])/( np.sqrt(combs["males_SE"]**2+combs["females_SE"]**2))) + math.log(2) ) / math.log(10))
else:
  combs["p_diff"]=0

combs.to_csv(prefix+"."+pheno+".sex_spec.gz", compression="gzip", sep=" ", index=False)

''' > source.py

# write checkpointing script
cat << "__EOF__" > checkpoint.sh
#!/bin/bash

while [ 1 -eq 1 ];do
    tar -cf checkpoint_new.tar *.regenie* *.log *.sex_spec.gz
    cp checkpoint.tar checkpoint_old.tar
    mv checkpoint_new.tar checkpoint.tar
    sleep 60
done
__EOF__
chmod +x checkpoint.sh

## check current progress
bgenix -g ${bgen} -list |grep -Ev "^#" > variants.list

mkdir checkpoint_folder
if test -f checkpoint.tar;then
    tar -xf checkpoint.tar -C checkpoint_folder/

fi
python3 checkpoint_status.py variants.list --analysis-type "${is_binary}" --endpoints "${sep="," phenolist}" --prefix "${prefix}"
# start checkpointing script
touch placeholder.regenie.placeholder placeholder.log placeholder.sex_spec.gz
./checkpoint.sh &
CHECKPOINT_PID=$!
## run regenie

# check if checkpoint file contains regenie outputs
if test -f phenos_done; then
    #copy files to here
    cp checkpoint_folder/*.regenie.gz ./
    cp checkpoint_folder/${prefix}.log ./
else 
    # test if we have processed some range
    if test -f remaining_range; then
        #move all non-complete male regenie outputs to have the range that they
        PROCESSED_RANGE=$(cat processed_range)
        for f in checkpoint_folder/*.regenie;do
            cp $f $f.$PROCESSED_RANGE
        done
        #copy all files with any range from checkpoint folder to main folder
        # this will fix the problem with one or more ranges missing when  
        cp checkpoint_folder/*.regenie.* ./
        cat checkpoint_folder/*.log* > ./${prefix}.log.$PROCESSED_RANGE
        REMAINING_RANGE=$(cat remaining_range)
        #we have, use the calculated remaining range 
        RANGE_OPTION="--range "$REMAINING_RANGE
        echo "Remaining range "$REMAINING_RANGE
    fi
    regenie \
    --step 2 \
    ${test_cmd} \
    ${if is_binary then "--bt --af-cc" else ""} \
    --bgen ${bgen} \
    --ref-first \
    --sample ${sample} \
    --covarFile ${cov_pheno} \
    --covarColList ${covariates} \
    --phenoFile ${cov_pheno} \
    --phenoColList ${sep="," phenolist} \
    --pred loco_list \
    ${if is_binary then "--use-null-firth firth_list" else ""} \
    --bsize ${bsize} \
    --threads $n_cpu \
    --out ${prefix} \
    $RANGE_OPTION \
    ${options}
    python3 combine_outputs.py "${prefix}*.regenie*" --output-name temp.regenie --endpoints "${sep="," phenolist}"
    #kill checkpointing for a moment. This is to prevent a situation where the partially compressed files are checkpointed, confusing the checkpointing script.
    kill $CHECKPOINT_PID 
    for p in ${sep=" " phenolist};do
        cat temp.regenie.$p |gzip >  ${prefix}"_"$p".regenie.gz"
    done
    #compress
    tar -cf checkpoint.tar *.regenie.gz ${prefix}.log*
    #restart checkpointing
    ./checkpoint.sh &
    CHECKPOINT_PID=$!
fi


## run sex-specific analysis

if [[ "${run_sex_specific}" == "true" ]];
then
    if test -f sex_spec_done; then
        #copy things 
        cp checkpoint_folder/*.sex_spec.gz ./
    else
        zcat ${cov_pheno} | awk 'NR==1{ for(i=1;i<=NF;i++) {h[$i]=i};
                                            if(!("${sex_col_name}" in h)) {
                                            print "Given sex column not found in phenotype file" > "/dev/stderr";
                                            exit 1;
                                            }
                                    }
                                NR>1&&$h["${sex_col_name}"]==0{ print $1,$1}' > males

        zcat ${cov_pheno} | awk 'NR==1{  for(i=1;i<=NF;i++) {h[$i]=i };
                                        if(!("${sex_col_name}" in h))
                                        {
                                        print "Given sex column not found in phenotype file" > "/dev/stderr";
                                        exit 1;
                                        }
                                    }
                                NR>1&&$h["${sex_col_name}"]==1{ print $1,$1}' > females

        sex_covars=$(echo ${covariates} | sed -e 's/${sex_col_name}//' | sed 's/^,//' | sed -e 's/,$//' | sed 's/,,/,/g')

        if [[ $sex_covars == ${covariates} ]];
        then
        echo "Warning! No sex covariate detected in used covariates."
        fi
        # checkpointing: Sex-spec phenos not done are in not_done_sex_spec,
        cat not_done_sex_spec | while read p || [[ -n $line ]];
        do
        #for p in ${sep=" " phenolist}; do

        if [[ "${is_binary}" == "true" ]];
        then
            N_cases_females=$(awk -v pheno=$p 'FNR==NR { females[$1]; next } FNR==1{ for(i=1;i<=NF;i++) {h[$i]=i} }; NR>1 && $h[pheno]==1 && $1 in females {print $1}' females <(zcat ${cov_pheno}) | wc -l)
            N_cases_males=$(awk -v pheno=$p 'FNR==NR { males[$1]; next } FNR==1{ for(i=1;i<=NF;i++) {h[$i]=i} }; NR>1 && $h[pheno]==1 && $1 in males {print $1}' males <(zcat ${cov_pheno}) | wc -l)
        else 
            N_cases_females=$(awk -v pheno=$p 'FNR==NR { females[$1]; next } FNR==1{ for(i=1;i<=NF;i++) {h[$i]=i} }; NR>1 && $h[pheno]!="NA" && $1 in females {print $1}' females <(zcat ${cov_pheno}) | wc -l)
            N_cases_males=$(awk -v pheno=$p 'FNR==NR { males[$1]; next } FNR==1{ for(i=1;i<=NF;i++) {h[$i]=i} }; NR>1 && $h[pheno]!="NA" && $1 in males {print $1}' males <(zcat ${cov_pheno}) | wc -l)
        fi 

        echo "Female cases: "$N_cases_females
        echo "Male cases: "$N_cases_males

        if [[ $N_cases_females -lt 10 ]] || [[ $N_cases_males -lt 10 ]];
        then
            echo "Less than 10 cases in a sex. Skipping testing of sex specific effects."
            touch ${prefix}"NOT_DONE.sex_spec.gz"
            continue
        fi

        base=${prefix}"_"$p".regenie.gz"

        zcat $base | awk 'NR==1{ for(i=1;i<=NF;i++) { h[$i]=i };
                                    if(!("ID" in h )||!("LOG10P" in h))
                                    {
                                    print "ID and LOG10P columns expected in association file" >"/dev/stderr"; exit 1}
                                    }
                            NR>1&&$h["LOG10P"]>${sex_specific_logpval} { print $h["ID"]} ' > $p".sex_variants"

        nvars=$(cat $p".sex_variants"| wc -l  )

        echo "$nvars variants will be tested for sex specific effects"

        if [[ $nvars -lt 1 ]];
        then
            ## NOTE
            ## NOTE: Echoing annoyingly the expected header here so that in gather the header can be taken from the first shard.
            ## NOTE This must match whats written from python below
        
            if [[ "${is_binary}" == "true" ]];
            then
                    echo -n "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ A1FREQ_CASES A1FREQ_CONTROLS INFO N TEST BETA SE CHISQ LOG10P"\
                    " EXTRA males_ID males_A1FREQ_CASES males_A1FREQ_CONTROLS males_N males_BETA males_SE males_LOG10P"\
                    " females_ID females_A1FREQ_CASES females_A1FREQ_CONTROLS females_N females_BETA females_SE females_LOG10P"\
                    " diff_beta p_diff" | bgzip > ${prefix}"."$p".sex_spec.gz"
                continue
            else 

                    echo -n "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P"\
                    " EXTRA males_ID males_N males_BETA males_SE males_LOG10P"\
                    " females_ID females_N females_BETA females_SE females_LOG10P"\
                    " diff_beta p_diff" | bgzip > ${prefix}"."$p".sex_spec.gz"
                continue
            fi 
        fi

        for s in males females;
        do
            echo "running $s analysis"
            echo "$(wc -l $s) individuals in file"
            head $s

            regenie \
            --step 2 \
            ${test_cmd} \
            ${if is_binary then "--bt --af-cc" else ""} \
            --bgen ${bgen} \
            --ref-first \
            --sample ${sample} \
            --keep $s  \
            --extract $p".sex_variants" \
            --covarFile ${cov_pheno} \
            --covarColList $sex_covars \
            --phenoFile ${cov_pheno} \
            --phenoColList $p \
            --pred ${pred} \
            --bsize ${bsize} \
            --threads $n_cpu \
            --gz \
            --out ${prefix}".sex_spec."$s \
            ${options}
        done

        echo "all files"
        ls *

        for f in $(ls *".sex_spec."*.regenie.gz ); do
            to=${DOLLAR}{f/%.regenie.gz/.gz}
            mv $f $to
            echo "file" $to
        done


        python3 source.py $p ${prefix} ${prefix}"_"$p".regenie.gz"
        # Checkpointing: recreate the checkpoint
        tar -cf checkpoint.tar *.regenie.gz ${prefix}.log* *.sex_spec.gz
        done
    fi
else
    ## sex specific analyses disabled but create the file so gather step is straightforward to implement
    touch ${prefix}"NOT_DONE.sex_spec.gz"
fi
rm placeholder.regenie.placeholder placeholder.log placeholder.sex_spec.gz
kill $CHECKPOINT_PID
    >>>

    output {
        Array[File] log = glob("*.log")
        Array[File] regenie = glob("*.regenie.gz")
        Array[File] sex_specifix = glob("*.sex_spec.gz")
    }

    runtime {
        docker: "${docker}"
        cpu: if length(phenolist) == 1 then 1 else if length(phenolist) <=4 then 2 else if length(phenolist) <= 10 then 4 else if length(phenolist) < 16 then 8 else 16
        memory:  "6 GB"
        disks: "local-disk " + (ceil(size(bgen, "G")) + 5) + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
        checkpointFile: "checkpoint.tar"
    }
}

task gather {

    Boolean is_binary
    Array[File] files

    Array[File] sex_files

    String docker

    command <<<

        set -euxo pipefail

        pheno=`basename ${files[0]} .regenie.gz | awk -F "." '{sub(/[^_]*_/, "", $NF); print $NF}'`
        mkdir regenie munged

        echo -e "`date`\tconcatenating result pieces into regenie/$pheno.regenie.gz, sorting by chr pos just in case"
        cat \
        <(zcat ${files[0]} | head -1) \
        <(for file in ${sep=" " files}; do
            zcat $file | tail -n+2
        done | sort -k1,1g -k2,2g) | bgzip > regenie/$pheno.regenie.gz


        if [[ ${length(sex_files)} -ge 1 ]]; then
          echo "concatenating sex diff results into regenie/$pheno.regenie.sex_diff.gz"
        cat \
        <(zcat ${sex_files[0]} | awk '{ print "#"$0; exit 0}') \
        <(for file in ${sep=" " sex_files}; do
            zcat $file | tail -n+2
        done | sort -k1,1g -k2,2g) | tr ' ' '\t' | bgzip > regenie/$pheno.regenie.sex_diff.gz
        tabix -s 1 -b 2 -e 2 regenie/$pheno.regenie.sex_diff.gz
        else
          touch regenie/$pheno.regenie.sex_diff.gz
          touch regenie/$pheno.regenie.sex_diff.gz.tbi
        fi

        if [[ "${is_binary}" == "true" ]]; then
            echo -e "`date`\tconverting to munged/$pheno.gz to a format used for importing to pheweb, omitting variants with -log10p NA (unsuccessful Firth/SPA)"
            zcat regenie/$pheno.regenie.gz | awk '
            BEGIN {FS=" "; OFS="\t"; split("CHROM GENPOS ALLELE0 ALLELE1 LOG10P BETA SE A1FREQ A1FREQ_CASES A1FREQ_CONTROLS", REQUIRED_FIELDS)}
            NR==1 {for(i=1;i<=NF;i++) h[$i]=i;
                   for(i in REQUIRED_FIELDS) if (!(REQUIRED_FIELDS[i] in h)) {print REQUIRED_FIELDS[i]" expected in regenie header">>"/dev/stderr"; exit 1}
                   print "#chrom","pos","ref","alt","pval","mlogp","beta","sebeta","af_alt","af_alt_cases","af_alt_controls"}
            NR>1 && $h["LOG10P"]!="NA" {print $h["CHROM"],$h["GENPOS"],$h["ALLELE0"],$h["ALLELE1"],10^-$h["LOG10P"],$h["LOG10P"],$h["BETA"],$h["SE"],$h["A1FREQ"],$h["A1FREQ_CASES"],$h["A1FREQ_CONTROLS"]}' \
            | bgzip > munged/$pheno.gz
        else
            echo -e "`date`\tconverting to munged/$pheno.gz to a format used for importing to pheweb"
            zcat regenie/$pheno.regenie.gz | awk '
            BEGIN {FS=" "; OFS="\t"; split("CHROM GENPOS ALLELE0 ALLELE1 LOG10P BETA SE A1FREQ", REQUIRED_FIELDS)}
            NR==1 {for(i=1;i<=NF;i++) h[$i]=i;
                   for(i in REQUIRED_FIELDS) if (!(REQUIRED_FIELDS[i] in h)) {print REQUIRED_FIELDS[i]" expected in regenie header">>"/dev/stderr"; exit 1}
                   print "#chrom","pos","ref","alt","pval","mlogp","beta","sebeta","af_alt"}
            NR>1  {print $h["CHROM"],$h["GENPOS"],$h["ALLELE0"],$h["ALLELE1"],10^-$h["LOG10P"],$h["LOG10P"],$h["BETA"],$h["SE"],$h["A1FREQ"]}' \
            | bgzip > munged/$pheno.gz
        fi

        echo -e "`date`\tprinting only chr, pos, pval to speed up and reduce memory use of qqplot.R"
        ## added SNP as updated manhattan function from qqman package fails without it.
        zcat munged/$pheno.gz | cut -f1,2,5 | awk 'BEGIN{ FS=OFS="\t"} NR==1{ print $0,"SNP"} NR>1{ print $0,NR}' > $pheno

        head $pheno

        echo -e "`date`\trunning qqplot.R"
        qqplot.R --file $pheno --chrcol "#chrom" --bp_col "pos" --pval_col "pval" --loglog_pval 10 > qq_out 2> qq_err

        echo -e "`date`\ttabixing munged/$pheno.gz"
        tabix -s1 -b2 -e2 munged/$pheno.gz
        echo -e "`date`\tdone"

    >>>

    output {
        File regenie = glob("regenie/*.regenie.gz")[0]
        File regenie_sex_diff = glob("regenie/*.regenie.sex_diff.gz")[0]
        File regenie_sex_diff_tbi = glob("regenie/*.regenie.sex_diff.gz.tbi")[0]
        File pheweb = glob("munged/*.gz")[0]
        File pheweb_tbi = glob("munged/*.gz.tbi")[0]
        File qq_out = "qq_out"
        File qq_err = "qq_err"
        Array[File] pngs = glob("*.png")
        Array[File] quantiles = glob("*qquantiles.txt")
    }

    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "8 GB"
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}

task summary{

    File input_file
    File finngen_annotation
    File finngen_tbi = finngen_annotation + ".tbi"

    Float summary_pval_thresh
    Float coding_pval_thresh
    String docker

    command <<<
        set -euxo pipefail
        python3 <<EOF
        #read file
        fname = "${input_file}"
        finngen_annotation_file = "${finngen_annotation}"
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
        FGAnnotation = NamedTuple('FGAnnotation',[('gene',str),('consequence',str),('rsid',str),('EXOME_enrichment_nfsee',str),('GENOME_enrichment_nfee',str)])
        def get_header(reader,file):
            with reader(file, "rt") as f:
                l = f.readline()
                return l.strip("\n").split("\t")

        def get_fg_annotation(iterator,variant,gene_idx,consequence_idx,rsid_idx,exome_enr_idx,genome_enr_idx,variant_idx) -> FGAnnotation:
            for v in iterator:
                cols = v.strip("\n").split("\t")
                if cols[variant_idx] == variant:
                    return FGAnnotation(cols[gene_idx],cols[consequence_idx],cols[rsid_idx],cols[exome_enr_idx],cols[genome_enr_idx])
            return FGAnnotation("","")

        #required columns
        fg_req_cols=["#variant","gene_most_severe","most_severe","rsid","EXOME_enrichment_nfsee","GENOME_enrichment_nfee"]
        #open finngen annotation tabix
        fg_tabix = pysam.TabixFile(finngen_annotation_file,parser=None)
        #get fg header column positions
        fg_header = get_header(gzip.open, finngen_annotation_file)
        fg_idx = {v:i for (i,v) in enumerate(fg_header)}
        #check for fg column existence
        if not all([a in fg_header for a in fg_req_cols]):
            raise Exception("Not all columns present in FinnGen annotation! Aborting...")
        var_idx, gene_idx, cons_idx, rsid_idx, exome_enr_idx, genome_enr_idx = (fg_idx["#variant"],fg_idx["gene_most_severe"],fg_idx["most_severe"],fg_idx["rsid"],fg_idx["EXOME_enrichment_nfsee"],fg_idx["GENOME_enrichment_nfee"])

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

                #add rsid, gene name, consequence, enrichments, phenotype
                header.extend(["rsid","gene_most_severe","most_severe","EXOME_enrichment_nfsee","GENOME_enrichment_nfee","phenotype"])
                summary_outfile.write("\t".join(header)+"\n")
                coding_outfile.write("\t".join(header)+"\n")

                #read lines
                for line in file:
                    line_columns = line.strip("\n").split('\t')
                    pvalue = float(line_columns[pval_idx])
                    if pvalue < coding_threshold or pvalue < summary_threshold:
                        cpra= (line_columns[cid],int(float(line_columns[pid])),line_columns[rid],line_columns[aid])
                        variant = "{}:{}:{}:{}".format(cpra[0],cpra[1],cpra[2],cpra[3])
                        fg_c = cpra[0].replace("chr","").replace("X","23").replace("Y","24").replace("XY","25").replace("MT","26").replace("M","26")
                        #annotate
                        fg_iter = fg_tabix.fetch(fg_c,cpra[1]-1, cpra[1])
                        fg_a = get_fg_annotation(fg_iter,variant,gene_idx,cons_idx,rsid_idx,exome_enr_idx,genome_enr_idx,var_idx)

                        line_columns.extend([
                            fg_a.rsid,
                            fg_a.gene,
                            fg_a.consequence,
                            fg_a.EXOME_enrichment_nfsee,
                            fg_a.GENOME_enrichment_nfee,
                            pheno,
                        ])
                        #gather row
                        #write to file
                        if pvalue < summary_threshold:
                            summary_outfile.write("\t".join(line_columns)+"\n")

                        #coding out
                        if pvalue < coding_threshold and fg_a.consequence in coding_groups:
                            coding_outfile.write("\t".join(line_columns)+"\n")
        print("summary created")
        EOF

    >>>

    output {
        File summary_out = glob("*_summary.txt")[0]
        File coding_out = glob("*_coding.txt")[0]
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
    Array[String] phenolist
    Boolean is_binary
    String cov_pheno
    String covariates

    String pred
    Array[String] loco
    String firth_list
    Array[String] nulls

    File bgenlist
    Array[String] bgens = read_lines(bgenlist)

    String sex_col_name

    Boolean is_single_sex
    Boolean run_sex_specific

    scatter (bgen in bgens) {
        call step2 {
            input: docker=docker, phenolist=phenolist, is_binary=is_binary, cov_pheno=cov_pheno,
              covariates=covariates, pred=pred, loco=loco, nulls=nulls,firth_list=firth_list,bgen=bgen,
              sex_col_name=sex_col_name, run_sex_specific=(run_sex_specific && !is_single_sex)
        }
    }

    Array[Array[String]] pheno_results = transpose(step2.regenie)

    Array[Array[String]] pheno_results_sex = transpose(step2.sex_specifix)

    scatter (pheno_result_idx in range(length(pheno_results))) {
        call gather {
            input: is_binary=is_binary, files=pheno_results[pheno_result_idx],
              sex_files = pheno_results_sex[pheno_result_idx], docker=docker
        }
        call summary {
            input: input_file=gather.pheweb, docker=docker
        }
    }

    output {
        Array[File] pheweb = gather.pheweb
        Array[File] pheweb_tbi = gather.pheweb_tbi
        Array[File] regenie_sex_diff = gather.regenie_sex_diff
        Array[File] regenie_sex_diff_tbi = gather.regenie_sex_diff_tbi
        Array[File] qq_out = gather.qq_out
        Array[File] qq_err = gather.qq_err
        Array[Array[File]] pngs = gather.pngs
        Array[Array[File]] quantiles = gather.quantiles
        Array[File] summary = summary.summary_out
        Array[File] coding = summary.coding_out
    }
}
