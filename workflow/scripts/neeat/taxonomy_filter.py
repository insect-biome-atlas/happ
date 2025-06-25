import pandas as pd
from argparse import ArgumentParser
import sys


def main(args):
    if not args.outfile:
        outfile = sys.stdout
    else:
        outfile = args.outfile
    taxdf = pd.read_csv(args.taxfile, sep="\t", index_col=0)
    taxdf = taxdf.loc[
        (~taxdf[args.rank].str.contains("_X+$"))
        & (~taxdf[args.rank].str.startswith("unclassified"))
    ]
    taxdf.to_csv(outfile, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-t",
        "--taxfile",
        type=str,
        help="Taxonomic and cluster assignment file",
        required=True,
    )
    parser.add_argument(
        "-r",
        "--rank",
        type=str,
        help="Rank for which to filter sequences",
        required=True,
    )
    parser.add_argument(
        "-o", "--outfile", type=str, help="Write filtered taxfile to outfile"
    )
    args = parser.parse_args()
    main(args)
