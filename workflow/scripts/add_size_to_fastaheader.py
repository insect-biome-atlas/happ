#!/usr/bin/env python
###################################

import pandas as pd
from Bio import SeqIO
from argparse import ArgumentParser
import sys


def read_sums(f):
    df = pd.read_csv(f, sep="\t", usecols=[0,1], names=["Seq","Sum"], header=0, index_col=0)
    return df.sort_values("Sum", ascending=False)


def read_seqs(f):
    seqs = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
    return seqs


def main(args):
    sums = read_sums(args.countsfile)
    seqs = read_seqs(args.seqsfile)
    with sys.stdout as fhout:
        for seqid in list(sums.index):
            c = sums.loc[seqid, "Sum"]
            try:
                fhout.write(f">{seqid};size={c}\n{str(seqs[seqid].seq)}\n")
            except KeyError:
                continue


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("seqsfile", type=str,
            help="Fasta file with sequences")
    parser.add_argument("countsfile", type=str,
            help="Tab-separated file with sums for sequences")
    args = parser.parse_args()
    main(args)
