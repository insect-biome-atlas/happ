#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
import numpy as np
import sys
from Bio.SeqIO import parse
import tqdm


def calc_counts(counts, method):
    if method == "sum":
        func = np.sum
    elif method == "median":
        func = np.median
    elif method == "mean":
        func = np.mean
    asv_counts = counts.groupby(level=0).sum().apply(func, axis=1)
    asv_counts = pd.DataFrame(asv_counts)
    asv_counts.columns = [method]
    return asv_counts


def get_rep(df, method):
    return df.loc[df[method] == df.max()[method]]


def filter_taxa(taxa, rank):
    filtered = taxa.loc[
        ~(taxa[rank].str.startswith("unclassified"))
        & ~(taxa[rank].str.contains("_[X]+$"))
    ]
    return filtered


def get_seqs(seqsfile, reps, ranks, rank):
    seqs = {}
    for record in parse(seqsfile, "fasta"):
        try:
            row = reps.loc[record.id]
        except KeyError:
            continue
        rank_name = row[rank]
        lineage = ";".join([x for x in row.loc[ranks]])
        desc = f"{lineage}"
        header = f">{record.id} {desc}"
        if rank_name not in seqs.keys():
            seqs[rank_name] = {"recid": header, "seq": record.seq}
        else:
            if len(record.seq) > len(seqs[rank_name]["seq"]):
                sys.stderr.write(
                    f"Replacing representative for {rank_name} with {record.id}\n"
                )
                seqs[rank_name] = {"recid": header, "seq": record.seq}
    return seqs


def write_seqs(seqs, outfile):
    if not outfile:
        w = sys.stdout
    else:
        w = open(outfile, "w")
    with w as fhout:
        for seqid, d in seqs.items():
            fhout.write(f"{d['recid']}\n{d['seq']}\n")


def main(args):
    taxa = pd.read_csv(args.taxa, index_col=0, sep="\t")
    taxa.index.name = "ASV"
    ranks = list(taxa.columns)
    filtered = filter_taxa(taxa, args.rank)
    counts = pd.read_csv(args.counts, index_col=0, sep="\t")
    counts.index.name = "ASV"
    asv_counts = calc_counts(counts, method=args.method)
    dataframe = pd.merge(filtered, asv_counts, left_index=True, right_index=True)
    reps = dataframe.groupby(args.rank).apply(get_rep, method=args.method)
    reps = reps.drop(args.rank, axis=1).reset_index()
    rep_size = reps.groupby(args.rank).size()
    sys.stderr.write(
        f"{rep_size.loc[rep_size>1].shape[0]} {args.rank} reps with >1 ASV\n"
    )
    reps.set_index("ASV", inplace=True)
    seqs = get_seqs(args.seqs, reps, ranks, args.rank)
    write_seqs(seqs, args.outfile)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("taxa", type=str, help="Taxonomic assignments")
    parser.add_argument(
        "counts", type=str, help="Counts of ASVs (rows) in samples (columns)"
    )
    parser.add_argument("seqs", type=str, help="Sequence file for ASVs")
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        help="Write representatives to outfile. Defaults to stdout",
    )
    parser.add_argument(
        "--rank", type=str, default="BOLD_bin", help="What level o group TAX ids"
    )
    parser.add_argument(
        "--method",
        type=str,
        choices=["sum", "median", "mean"],
        default="sum",
        help="Method to select representative base on counts. Defaults to 'sum' across samples",
    )
    args = parser.parse_args()
    main(args)
