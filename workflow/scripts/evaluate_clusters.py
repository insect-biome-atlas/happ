#!/usr/bin/env python

import pandas as pd
from argparse import ArgumentParser
import sys
from tqdm import tqdm
from sklearn.metrics import homogeneity_score as homscore
from sklearn.metrics import completeness_score as comscore


def read_asv_clusters(files):
    """
    Read cluster files
    """
    res = pd.DataFrame()
    for f in files:
        _ = pd.read_csv(f, sep="\t", index_col=0, dtype=str)
        # Only take first ASV/cluster combination. There shouldn't be more
        # than one but just in case
        _ = _.groupby(level=0).first()
        res = pd.concat([_, res])
    return res


def merge_cluster_df(clustdf, taxdf):
    """
    Merge cluster dataframe with taxonomy
    """
    dataf = pd.merge(clustdf, taxdf, left_index=True, right_index=True, how="inner")
    dataf["cluster_Family"] = dataf["cluster"] + dataf["Family"]
    return dataf


def pairs(df):
    """
    Return all possible pairs from given dataframe
    """
    N = float(df.shape[0])
    return N * (N - 1) / 2


def truepos(res, rank):
    """
    Take dataframe as input and groupby rank, then return number of pairs
    """
    return res.groupby(rank).apply(pairs)


def falseNegatives(df, cluster_col, rank, silent=False):
    """
    Iterates each unique label (e.g. species) and calculates how many that should
    be clustered but are not
    """
    cl_rank_size = pd.DataFrame(df.groupby([cluster_col, rank]).size()).reset_index()
    FN = 0
    if not silent:
        iter = tqdm(
            cl_rank_size[rank].unique(),
            desc="finding false negatives",
            ncols=10,
            dynamic_ncols=True,
            unit=f" {rank}",
        )
    else:
        iter = cl_rank_size[rank].unique()
    for doc in iter:
        items = list(cl_rank_size.loc[cl_rank_size[rank] == doc, 0].items())
        for i, item in enumerate(items):
            FN += item[1] * sum([x[1] for x in items[i + 1 :]])
    return FN


def precision_recall(df, cluster_col, rank, silent=False):
    """
    Calculates true positives, false positives, false negatives and returns
    precision and recall
    """
    totalClusters = len(df[cluster_col].unique())
    totalTaxa = len(df[rank].unique())
    totalPositives = sum(df.groupby(cluster_col).apply(pairs))
    if totalTaxa == 1:
        TP = df.groupby(cluster_col).apply(truepos, rank=rank).sum().values[0]
    else:
        TP = df.groupby(cluster_col).apply(truepos, rank=rank).sum()
    FP = totalPositives - TP
    # totalNegatives = pairs(df) - totalPositives
    FN = falseNegatives(df, cluster_col, rank, silent)
    # TN = totalNegatives - FN
    if TP + FP == 0:
        precision = 0.0
    else:
        precision = float(TP) / (TP + FP)
    if TP + FN == 0:
        recall = 0.0
    else:
        recall = float(TP) / (TP + FN)
    if not silent:
        sys.stdout.write(f"Total clusters {totalClusters}\n")
        sys.stdout.write(f"Total number of {rank} {totalTaxa}\n")
        sys.stdout.write(f"Total positives: {totalPositives}\n")
        sys.stdout.write(f"True positives: {TP}\n")
        sys.stdout.write(f"False positives: {FP}\n")
        sys.stdout.write(f"False negatives: {FN}\n")
    return precision, recall


def homcom(clustdf, clust_col, rank):
    labels_true = clustdf[rank]
    labels_pred = clustdf[clust_col]
    homogeneity = homscore(labels_true, labels_pred)
    completeness = comscore(labels_true, labels_pred)
    return homogeneity, completeness


def calc_order_level(clustdf, rank):
    order_level = {}
    orders = clustdf.Order.unique()
    for o in tqdm(
        orders,
        desc="calculating stats",
        ncols=10,
        dynamic_ncols=True,
        unit=f" Order",
    ):
        _dataf = clustdf.loc[clustdf.Order == o]
        totalClusters = len(_dataf["cluster_Family"].unique())
        n_rank = len(_dataf[rank].unique())
        n_asv = _dataf.shape[0]
        p, r = precision_recall(_dataf, "cluster_Family", rank, silent=True)
        h, c = homcom(_dataf, "cluster_Family", rank)
        order_level[o] = {
            "precision": p,
            "recall": r,
            "homogeneity": h,
            "completeness": c,
            rank: n_rank,
            "ASVs": n_asv,
            "clusters": totalClusters,
        }
    order_leveldf = pd.DataFrame(order_level).T
    order_leveldf.index.name = "Order"
    return order_leveldf


def main(args):
    # Read taxonomies
    asv_taxa = pd.read_csv(args.taxfile, sep="\t", index_col=0)
    sys.stderr.write(f"#{asv_taxa.shape[0]} ASV taxonomies loaded\n")
    # Remove ASVs without assignments for args.rank
    sys.stderr.write(f"#Removing ASVs without assignments for {args.rank}\n")
    asv_taxa = asv_taxa.loc[
        (~asv_taxa[args.rank].str.contains("_X+$"))
        & (~asv_taxa[args.rank].str.startswith("unclassified"))
    ]
    sys.stderr.write(f"#{asv_taxa.shape[0]} ASVs remaining\n")
    # Read cluster files
    sys.stderr.write(f"#Loading cluster results from {len(args.clustfiles)} files\n")
    clustdf = read_asv_clusters(args.clustfiles)
    # Merge with taxonomies
    sys.stderr.write("#Merging with taxonomic assignments\n")
    clustdf = merge_cluster_df(clustdf, asv_taxa)
    sys.stderr.write(f"#{clustdf.shape[0]} ASVs remaining after merging\n")
    if clustdf.shape[0] == 0:
        sys.exit("Not enough data to evaluate\n")
    if args.order_level:
        sys.stderr.write(
            f"Calculating statistics for {len(clustdf.Order.unique())} orders\n"
        )
        order_leveldf = calc_order_level(clustdf, args.rank)
        sys.stderr.write(f"Writing order-level results to {args.order_level}\n")
        order_leveldf.to_csv(args.order_level, sep="\t")
    precision, recall = precision_recall(clustdf, "cluster_Family", args.rank)
    homogeneity, completeness = homcom(clustdf, "cluster_Family", args.rank)
    sys.stdout.write(f"precision\t{precision}\nrecall\t{recall}\n")
    sys.stdout.write(f"homogeneity\t{homogeneity}\ncompleteness\t{completeness}\n")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "taxfile", type=str, help="File with taxonomic assignments for ASVs"
    )
    parser.add_argument(
        "clustfiles",
        type=str,
        nargs="+",
        help="Cluster membership file. Should be tab-separated,"
        " have a header, and contain ASV ids in the first column and "
        "cluster name in second column",
    )
    parser.add_argument(
        "--rank",
        type=str,
        default="Species",
        help="Evaluate clusters against a 'true' cluster rank",
    )
    parser.add_argument(
        "--order_level",
        type=str,
        help="Calculate order-level stats and write to this file",
    )
    args = parser.parse_args()
    main(args)
