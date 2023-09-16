#!/usr/bin/env python

import pandas as pd
from argparse import ArgumentParser
import gzip as gz
import tqdm
import sys
from collections import defaultdict


def get_pdist(df, rank, directory):
    reps = df.loc[df.representative==1]
    pdist = {}
    rep_asvs = list(reps.index)
    rep_dict = defaultdict.fromkeys(rep_asvs)
    for tax in tqdm.tqdm(sorted(reps[rank].unique()), desc="reading distfiles", position=0, leave=False, ncols=100, unit=f"{rank}s"):
        distfile = f"{directory}/{tax}/asv_seqs.dist.gz"
        _pdist = {}
        with gz.open(distfile, 'rt') as fhin:
            for line in tqdm.tqdm(fhin, desc=tax, position=1, leave=False, ncols=100, unit=" lines"):
                asv1, asv2, pident = line.split("\t")
                pident = float(pident)
                try:
                    rep_dict[asv1]
                    rep_dict[asv2]
                except KeyError:
                    continue
                try:
                    _pdist[asv1][asv2] = pident
                except KeyError:
                    _pdist[asv1] = {asv2: pident}
        pdist.update(_pdist)
    return pdist


def main(args):
    df = pd.read_csv(args.clustfile, sep="\t", index_col=0, nrows=args.nrows)
    rank = args.rank
    directory = args.directory
    pdist = get_pdist(df, rank, directory)
    for key, d in pdist.items():
        for key2, val in d.items():
            sys.stdout.write(f"{key}\t{key2}\t{val}\n")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("clustfile", type=str,
                        help="File mapping asvs to clusters, and taxonomy")
    parser.add_argument("--directory", type=str,
                        help="Main directory under which vsearch distance output is stored")
    parser.add_argument("--rank", type=str, default="Family",
                        help="Rank at which data has been split")
    parser.add_argument("--nrows", type=int,
                        help="Number of rows to read from clustfile (debugging purposes only)")
    args = parser.parse_args()
    main(args)
