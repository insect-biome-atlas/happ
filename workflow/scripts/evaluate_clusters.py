#!/usr/bin/env python

import pandas as pd
from argparse import ArgumentParser
import sys


def main(args):
    taxdf = pd.read_csv(args.taxfile, sep="\t", header=0, index_col=0)
    clustdf = pd.read_csv(args.clustfile, sep="\t", header=0, index_col=0)
    merged = pd.merge(clustdf, taxdf, left_index=True, right_index=True)
    nunique = merged.groupby("cluster").nunique()
    nunique.to_csv(sys.stdout, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("taxfile", type=str,
                        help="File with taxonomic assignments for ASVs")
    parser.add_argument("clustfile", type=str,
                        help="Cluster membership file. Should be tab-separated,"
                             " have a header, and contain ASV ids in the first column and "
                             "cluster name in second column")
    args = parser.parse_args()
    main(args)