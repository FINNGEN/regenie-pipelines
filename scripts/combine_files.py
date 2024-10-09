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