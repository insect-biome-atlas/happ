#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
import numpy as np
import sys
import os
from Bio.SeqIO import parse
import tqdm
import time


def rank_max(df):
    idx = df.idxmax(numeric_only=True)
    for method in ["median", "mean", "sum"]:
        if len(df[method].unique()) > 1:
            break
    rep = idx.loc[method]
    return pd.DataFrame(df.loc[rep]).T


def get_reps(df, rank):
    start = time.time()
    reps = df.groupby(rank, group_keys=False).apply(rank_max)
    sys.stderr.write(f"Found representatives in {time.time()-start} seconds\n")
    return reps


def filter_taxa(taxa, rank):
    filtered = taxa.loc[
        ~(taxa[rank].str.startswith("unclassified"))
        & ~(taxa[rank].str.contains("_[X]+$"))
    ]
    return filtered


def get_seqs(seqsfile, reps, ranks, rank):
    start = time.time()
    seqs = {}
    for record in tqdm.tqdm(parse(seqsfile, "fasta"), unit=" records",):
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
                seqs[rank_name] = {"recid": header, "seq": record.seq}
    sys.stderr.write(f"Read seqs in {time.time()-start} seconds\n")
    return seqs


def write_seqs(seqs, outfile):
    if not outfile:
        w = sys.stdout
    else:
        w = open(outfile, "w")
    with w as fhout:
        for seqid, d in seqs.items():
            fhout.write(f"{d['recid']}\n{d['seq']}\n")


def calc_colsum(f):
    colsums = []
    data = pd.read_csv(f, sep="\t", index_col=0, chunksize=100000)
    for item in tqdm.tqdm(
        data, unit=" chunks", ncols=50, desc="Reading counts file in chunks"
    ):
        colsums.append(item.sum(axis=0))
    return sum(colsums)


def read_input(f):
    start = time.time()
    df = pd.read_csv(f, index_col=0, sep="\t")
    sys.stderr.write(f"Read {f} in {time.time()-start} seconds\n")
    return df


def calc_counts(countsdf, colsums, ids=None, normalize=False):
    if ids is None:
        ids = []
    if len(ids) > 0:
        common = list(set(countsdf.index).intersection(ids))
        countsdf = countsdf.loc[common]
    countsdf.fillna(0, inplace=True)
    start = time.time()
    if normalize:
        size_factor = colsums.loc[countsdf.columns][colsums.columns[0]]
        countsdf = countsdf.div(size_factor, axis=1)
    median = pd.DataFrame(countsdf.apply(np.nanmedian, axis=1), columns=["median"])
    mean = pd.DataFrame(countsdf.apply(np.mean, axis=1), columns=["mean"])
    s = pd.DataFrame(countsdf.apply(np.sum, axis=1), columns=["sum"])
    calculated_counts = pd.merge(
        pd.merge(median, mean, left_index=True, right_index=True),
        s, left_index=True, right_index=True
    )
    
    sys.stderr.write(
        f"Calculated counts in" f"{time.time()-start} seconds\n"
    )
    return calculated_counts


def main(args):
    sys.stderr.write(f"Reading taxfile {args.taxa}\n")
    taxa = read_input(args.taxa)
    if args.prefix:
        taxa[args.rank] = [f"{args.prefix}_{x}" for x in taxa[args.rank]]
    sys.stderr.write(f"Read {taxa.shape[0]} records\n")
    taxa.index.name = "ASV"
    ranks = list(taxa.columns)
    if args.no_unclassified:
        sys.stderr.write(f"Filtering taxonomy to remove unassigned\n")
        filtered = filter_taxa(taxa, args.rank)
        sys.stderr.write(f"{filtered.shape[0]} records remaining\n")
    else:
        filtered = taxa
    colsums = None
    if args.normalize:
        if args.colsums:
            sys.stderr.write(f"Reading column sums from {args.colsums}\n")
            colsums = read_input(args.colsums)
        else:
            sys.stderr.write(f"Calculating column sums for normalization\n")
            colsums = calc_colsum(args.counts)
    sys.stderr.write(
        f"Reading countsfile {args.counts} and calculating abundance of ASVs "
        f"using {args.method} across samples\n"
    )
    countsdf = read_input(args.counts)
    counts = calc_counts(
        countsdf,
        colsums,
        ids=list(filtered.index),
        normalize=args.normalize,
    )
    counts.index.name = "ASV"
    dataframe = pd.merge(filtered, counts, left_index=True, right_index=True)
    sys.stderr.write(f"Finding representatives for rank {args.rank}\n")
    if dataframe.shape[0] == 1:
        reps = dataframe
    else:
        reps = get_reps(dataframe, args.rank)
    rep_size = reps.groupby(args.rank).size()
    sys.stderr.write(
        f"{rep_size.loc[rep_size>1].shape[0]} {args.rank} reps with >1 ASV\n"
    )
    if args.taxa_table:
        sys.stderr.write(f"Adding taxonomic info from {args.taxa_table}\n")
        extra_taxdf = read_input(args.taxa_table)
        taxdf = pd.merge(dataframe, extra_taxdf, left_index=True,
                         right_index=True)
        taxdf = taxdf.assign(
            representative=pd.Series([0] * taxdf.shape[0], index=taxdf.index)
        )
        taxdf.loc[reps.index, "representative"] = 1
        taxaout = f"{os.path.dirname(args.outfile)}/cluster_taxonomy.tsv"
        taxdf.to_csv(taxaout, sep="\t")
    sys.stderr.write(f"Reading sequences from {args.seqs}\n")
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
        "outfile",
        type=str,
        help="Write representatives to outfile",
    )
    parser.add_argument(
        "--rank", type=str, default="BOLD_bin", help="What level to group TAX "
                                                     "ids"
    )
    parser.add_argument("--prefix", type=str, help="Prefix cluster name with "
                                                   "string")
    parser.add_argument(
        "--method",
        type=str,
        choices=["sum", "median", "mean"],
        default="median",
        help="Method to select representative base on counts. Defaults to "
             "'sum' across samples",
    )
    parser.add_argument(
        "--normalize",
        action="store_true",
        help="Normalize counts before applying method",
    )
    parser.add_argument(
        "--colsums",
        type=str,
        help="File with pre-calculated read count sums for samples",
    )
    parser.add_argument(
        "--no-unclassified",
        action="store_true",
        help="Remove ASVs marked as 'unclassified' or suffixed with '_X'",
    )
    parser.add_argument(
        "--taxa-table",
        type=str,
        help="Additional taxonomy table to merge representatives with",
    )
    args = parser.parse_args()
    main(args)
