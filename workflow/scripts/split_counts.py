#!/usr/bin/env python

import pandas as pd
from argparse import ArgumentParser
import polars as pl


def main(args):
    r = pd.read_csv(args.counts, sep="\t", nrows=1, index_col=0)
    df = pl.read_csv(args.counts, has_header=True, sep="\t")
    samples = list(r.columns)
    for sample in samples:
        outfile=f"{args.outdir}/{sample.replace(' ', '_')}.sum.tsv"
        _ = df.filter(pl.col(sample)>0).select(["ASV_ID", sample])
        _.columns = list(map(lambda x: x.replace(sample, "Sum"), _.columns))
        _.write_csv(outfile, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "counts", type=str, help="Counts file for sequences (rows) in samples (columns)"
    )
    parser.add_argument("outdir", type=str, help="Output directory")
    args = parser.parse_args()
    main(args)
