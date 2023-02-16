#!/usr/bin/env python

import pandas as pd
from argparse import ArgumentParser


def main(args):
    r = pd.read_csv(args.counts, sep="\t", nrows=1, index_col=0)
    samples = list(r.columns)
    for i, sample in enumerate(samples, start=1):
        df = pd.read_csv(args.counts, usecols=[0, i], sep="\t", index_col=0)
        df = df.loc[df[sample] > 0].rename(columns={sample: "Sum"})
        df.to_csv(f"{args.outdir}/{sample}.sum.tsv", sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "counts", type=str, help="Counts file for sequences (rows) in samples (columns)"
    )
    parser.add_argument("outdir", type=str, help="Output directory")
    args = parser.parse_args()
    main(args)
