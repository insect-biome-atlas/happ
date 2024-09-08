#!/usr/bin/env python

from argparse import ArgumentParser
from Bio.SeqIO import parse
import pandas as pd

def main(args):
    records = []
    with open(args.fasta, "r") as fhin, open(args.out_prefix + "_filtered.fna", "w") as fhout:
        for record in parse(fhin, "fasta"):
            if len(record.seq) >= args.min_length:
                if args.max_length is not None:
                    if len(record.seq) <= args.max_length:
                        records.append(record.id)
                        fhout.write(f">{record.id}\n{record.seq}\n")
    if args.count_table:
        with open(args.count_table, "r") as fhin, open(args.out_prefix + "_filtered.table.tsv", "w") as fhout:
            for i, line in enumerate(fhin):
                if i==0:
                    fhout.write(line)
                    continue
                asv = line.strip().split("\t")[0]
                if asv in records:
                    fhout.write(line)
            

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-f",
        "--fasta",
        help="Input fasta file",
        required=True,
    )
    parser.add_argument(
        "-m",
        "--min_length",
        help="Minimum length of sequence to keep",
        required=True,
        type=int,
        default=0,
    )
    parser.add_argument(
        "-M",
        "--max_length",
        help="Maximum length of sequence to keep",
        required=False,
        type=int,
    )
    parser.add_argument(
        "-t", 
        "--count_table",
        help="Count table file of the ASVs ",
    )
    parser.add_argument(
        "-p",
        "--out_prefix",
        help="Prefix of the output files with filtered ASVs and the corresponding count table. The output files will be in *<prefix>_filtered.* format",
        required=True,
    )
    args = parser.parse_args()
    main(args)