#!/usr/bin/env python

import pandas as pd
import os
import subprocess
from argparse import ArgumentParser
import sys
from tempfile import NamedTemporaryFile

def main(args):
    countsfile = args.countsfile
    sampleid = args.sampleid
    _ = pd.read_csv(countsfile, sep="\t", header=0, nrows=1)
    header = _.columns.tolist()
    sample_index = header.index(sampleid) + 1
    sys.stderr.write(f"Reading counts for sample {sampleid} from {countsfile}\n")
    with NamedTemporaryFile(delete=False) as fhout:
        subprocess.run(["cut",f"-f1,{sample_index}",countsfile], stdout=fhout)
        fhout.close()
        with open(fhout.name, 'r') as fhin:
            for i, line in enumerate(fhin):
                line = line.rstrip()
                if i==0:
                    line = line.replace(sampleid, "Sum")
                    print(line)
                    continue
                count = line.split("\t")[1]
                if int(count) > 0:
                    print(line)
    os.remove(fhout.name)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("countsfile", type=str, help="File with counts for ASVs (rows) in samples (columns)")
    parser.add_argument("sampleid", type=str, help="Sample id to extract counts for")
    args = parser.parse_args()
    main(args)
    
