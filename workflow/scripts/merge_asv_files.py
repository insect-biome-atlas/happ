#!/usr/bin/env python

from argparse import ArgumentParser
from Bio.SeqIO import parse
import pandas as pd
import sys
from tqdm import tqdm


def read_seqs(files):
    seqs = {}
    idmap = {}
    for f in files:
        sys.stderr.write(f"Reading seqs from {f}\n")
        for record in tqdm(parse(f, "fasta"), unit=" seqs", ncols=100):
            try:
                index = seqs[str(record.seq)]
            except KeyError:
                index = record.id
                seqs[str(record.seq)] = index
            idmap[record.id] = index
    return seqs, idmap


def read_counts(files, idmap):
    counts = pd.DataFrame()
    for f in files:
        sys.stderr.write(f"Reading counts from {f}\n")
        _counts = pd.read_csv(f, sep="\t", index_col=0)
        _counts = _counts.rename(index=idmap)
        counts = pd.merge(
            counts, _counts, left_index=True, right_index=True, how="outer"
        )
    return counts.fillna(0)


def main(args):
    seqs, idmap = read_seqs(args.fasta)
    counts = read_counts(args.count, idmap)
    sys.stderr.write(f"Writing merged counts to {args.counts_out}\n")
    counts.to_csv(args.counts_out, sep="\t")
    sys.stderr.write(f"Writing sequences to {args.fasta_out}\n")
    with open(args.fasta_out, "w") as fhout:
        for seq, seqid in seqs.items():
            fhout.write(f">{seqid}\n{seq}\n")
    sys.stderr.write(f"Writing seqid map to {args.idmap}\n")
    with open(args.idmap, "w") as fhout:
        for seq1, seq2 in idmap.items():
            fhout.write(f"{seq1}\t{seq2}\n")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--fasta", nargs="+", help="Fasta files")
    parser.add_argument("--count", nargs="+", help="Counts files")
    parser.add_argument("--fasta_out", type=str, help="Seq file output")
    parser.add_argument("--counts_out", type=str, help="Merged table output")
    parser.add_argument(
        "--idmap", type=str, help="Mapping table showing which ASV ids share sequence"
    )
    args = parser.parse_args()
    main(args)
