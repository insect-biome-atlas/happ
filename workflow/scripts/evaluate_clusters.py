#!/usr/bin/env python

import pandas as pd
from argparse import ArgumentParser
import sys


def read_clustfiles(files):
    clusters = {}
    clust_num = 0
    for file in files:
        with open(file, 'r') as fhin:
            file_clusters = {}
            for i, line in enumerate(fhin):
                if i == 0:
                    continue
                asv, cluster = line.rstrip().rsplit()
                if cluster not in file_clusters.keys():
                    clust_num += 1
                    cluster_name = f"cluster{clust_num}"
                    file_clusters[cluster] = cluster_name
                else:
                    cluster_name = file_clusters[cluster]
                clusters[asv] = cluster_name
    return pd.DataFrame(clusters, index=["cluster"]).T



def main(args):
    taxdf = pd.read_csv(args.taxfile, sep="\t", header=0, index_col=0)
    clustdf = read_clustfiles(args.clustfiles)
    merged = pd.merge(clustdf, taxdf, left_index=True, right_index=True)
    nunique = merged.groupby("cluster").nunique()
    nunique.to_csv(sys.stdout, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("taxfile", type=str,
                        help="File with taxonomic assignments for ASVs")
    parser.add_argument("clustfiles", type=str, nargs="+",
                        help="Cluster membership file. Should be tab-separated,"
                             " have a header, and contain ASV ids in the first column and "
                             "cluster name in second column")
    args = parser.parse_args()
    main(args)