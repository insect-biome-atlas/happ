#!/usr/bin/env python
###################################

import polars as pl
import pandas as pd
from argparse import ArgumentParser
import sys


def main(args):
    r = pd.read_csv(args.infile, sep="\t", nrows=1, index_col=0)
    index_name = r.index.name
    df = pl.read_csv(args.infile, has_header=True, separator="\t")
    asv_tot = df.select(pl.col(pl.Int64)).sum(axis=1)
    x = df.select(index_name).to_series()
    asvs = [x[i] for i in range(0, x.shape[0])]
    asv_summed = pd.DataFrame(
        {"ASV_ID": asvs, "Sum": [asv_tot[i] for i in range(0, asv_tot.shape[0])]}
    )
    asv_summed.set_index("ASV_ID", inplace=True)
    asv_summed.to_csv(sys.stdout, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("infile", type=str, help="Input file with counts for ASVs")
    args = parser.parse_args()
    main(args)
