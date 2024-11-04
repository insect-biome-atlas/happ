#!/usr/bin/env python

from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re


def read_records(infile, input_format, truncate=False, mapfile=None):
    m = {}
    records = []
    with open(infile, 'r') as fhin:
        for i, record in enumerate(SeqIO.parse(fhin, input_format), 1):
            m[record.id] = f"rec{i}"
            if truncate:
                record.id = f"rec{i}"
            records.append(record)

    if truncate:
        with open(mapfile, 'w') as fhout:
            for key, value in m.items():
                fhout.write(f"{key}\t{value}\n")
    return records


def main(args):
    mapfile = f"{args.outfile}.map"
    truncate = False
    if args.output_format == "phylip":
        truncate = True
    records = read_records(infile=args.infile, input_format=args.input_format, truncate=truncate, mapfile=mapfile)

    count = SeqIO.write(records, args.outfile, args.output_format)
    print(f"Converted {count} records")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("infile", type=str,
                        help="Input file alignment ")
    parser.add_argument("input_format", type=str,
                        choices=["nexus","stockholm","fasta","phylip"],
                        help="Format of input alignment")
    parser.add_argument("outfile", type=str,
                        help="Output alignment")
    parser.add_argument("output_format", type=str,
                        choices=["nexus", "stockholm", "fasta", "phylip"],
                        help="Format of output alignment")
    args = parser.parse_args()
    main(args)