#!/usr/bin/env python

import pandas as pd
from argparse import ArgumentParser

def main(args):
    taxonomy = pd.read_csv(args.taxonomy, sep="\t", index_col=0)
    taxdf = taxonomy.loc[taxonomy[args.rank]==args.taxa]
    counts = pd.read_csv(args.counts, sep="\t", index_col=0)
    if "representative" in taxdf.columns:
        taxdf = taxdf.loc[taxdf.representative==1]
    else:
        taxdf["representative"] = [1]*taxdf.shape[0]
    if "cluster" in taxdf.columns:
        seqs = list(set(taxdf["cluster"]).intersection(counts.index))
        taxdf = taxdf.set_index("cluster")
    else:
        seqs = list(set(taxdf.index).intersection(counts.index))
    tax_counts = counts.loc[seqs]
    tax_counts.to_csv(f"{args.outdir}/{args.taxa}_counts.tsv", sep="\t")
    taxdf.to_csv(f"{args.outdir}/{args.taxa}_taxonomy.tsv", sep="\t")
    

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("taxonomy", type=str, help="Taxonomy file")
    parser.add_argument("counts", type=str, help="Counts file")
    parser.add_argument("-r", "--rank", type=str, help="Taxonomic rank", required=True)
    parser.add_argument("-t", "--taxa", type=str, help="Taxa name", required=True)
    parser.add_argument("-o", "--outdir", type=str, help="Output directory", required=True)
    args = parser.parse_args()
    main(args)
