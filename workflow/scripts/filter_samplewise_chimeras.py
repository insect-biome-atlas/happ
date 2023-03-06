#!/usr/bin/env python

from argparse import ArgumentParser
from Bio.SeqIO import parse
import sys
import gzip as gz
import numpy as np


def main(args):
    chims = {}
    sys.stderr.write(f"Reading chimeras from {len(args.chims)} files\n")
    for f in args.chims:
        with gz.open(f, "rt") as fhin:
            for record in parse(fhin, "fasta"):
                seqid, _ = (record.id).split(";")
                size = _.split("=")[1]
                try:
                    chims[seqid].append(int(size))
                except KeyError:
                    chims[seqid] = [int(size)]
    sys.stderr.write(f"Read {len(chims.keys())} chimeras\n")
    sys.stderr.write(f"Reading all sequences from {args.fasta}\n")
    nonchims = 0
    with sys.stdout as fhout:
        for record in parse(args.fasta, "fasta"):
            try:
                chims[record.id]
            except KeyError:
                fhout.write(f">{record.id}\n{record.seq}\n")
                nonchims += 1
    sys.stderr.write(f"Wrote {nonchims} non-chimeric seqs\n")
    if args.chimeraids:
        with open(args.chimeraids, "w") as fhout:
            fhout.write("ASV\tn_samples\tmean_size\tmedian_size\n")
            for seqid, l in chims.items():
                fhout.write(f"{seqid}\t{len(l)}\t{np.mean(l)}\t{np.median(l)}\n")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--chims", nargs="+", help="Fasta file(s) with chimeras")
    parser.add_argument("--fasta", type=str, help="Fasta file of all seqs")
    parser.add_argument("--chimeraids", type=str, help="Write chimeric ids to file")
    args = parser.parse_args()
    main(args)
