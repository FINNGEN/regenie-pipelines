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