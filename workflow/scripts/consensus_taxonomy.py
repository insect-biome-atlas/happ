#!/usr/bin/env python
import argparse
from argparse import ArgumentParser
import pandas as pd
import tqdm
import sys
from collections import defaultdict


def find_consensus_taxonomies(
    clustdf, clust_column, ranks, consensus_ranks, consensus_threshold
):
    cluster_taxonomies = {}
    cons_ranks_reversed = consensus_ranks.copy()
    cons_ranks_reversed.reverse()
    for cluster in tqdm.tqdm(
        clustdf[clust_column].unique(),
        desc="finding consensus taxonomies",
        unit=" clusters",
    ):
        rows = clustdf.loc[clustdf[clust_column] == cluster]
        lineage = {cluster: defaultdict(lambda: "unresolved")}
        lineage[cluster].update({rank: "unresolved" for rank in ranks})
        taxlabel = "unresolved"
        for rank in cons_ranks_reversed:
            # Sum ASV sums up to rank
            rank_sums = rows.groupby(rank).sum(numeric_only=True)
            # Calculate percent for rank labels
            rank_sums_percent = rank_sums.div(rank_sums.sum()) * 100
            # Sort percent values in descending order
            rank_sums_percent.sort_values("ASV_sum", ascending=False, inplace=True)
            # Create a list <above_thresh> of number of labels at or above threshold
            above_thresh = list(
                rank_sums_percent.loc[
                    rank_sums_percent["ASV_sum"] >= consensus_threshold
                ].index
            )
            # If only one assignment is above threshold, use this lineage to resolve taxonomy
            if len(above_thresh) == 1:
                taxlabel = above_thresh[0]
                lineage = (
                    rows.loc[rows[rank] == above_thresh[0]][ranks]
                    .head(1)
                    .to_dict(orient="index")
                )
                break
        cluster_taxonomies[cluster] = list(lineage.values())[0]
        ranks_below = cons_ranks_reversed[0 : cons_ranks_reversed.index(rank)]
        if taxlabel == "unresolved":
            prefix = ""
        else:
            prefix = "unresolved."
        for r in ranks_below:
            cluster_taxonomies[cluster][r] = f"{prefix}{taxlabel}"
    return pd.DataFrame(cluster_taxonomies).T


def main(args):
    sys.stderr.write("####\n" f"Reading ASV clusters from {args.clustfile}\n")
    clustdf = pd.read_csv(args.clustfile, sep="\t", index_col=0, header=0)
    sys.stderr.write("####\nSumming counts for ASVs\n")
    asv_sum = pd.read_csv(args.countsfile, sep="\t", index_col=0, header=None, names=["ASV", "ASV_sum"])
    clustdf = clustdf.loc[:, [args.clust_column] + args.ranks]
    clustdf = pd.merge(asv_sum, clustdf, left_index=True, right_index=True)
    sys.stderr.write(
        "####\n"
        f"Read {clustdf.shape[0]} ASVs in {len(clustdf[args.clust_column].unique())} clusters\n"
    )
    sys.stderr.write(
        "####\n"
        f"Resolving taxonomies using {args.consensus_threshold}% majority rule threshold\n"
    )
    resolved = find_consensus_taxonomies(
        clustdf=clustdf,
        clust_column=args.clust_column,
        ranks=args.ranks,
        consensus_ranks=args.consensus_ranks,
        consensus_threshold=args.consensus_threshold,
    )
    resolved.index.name = "cluster"
    resolved.sort_index(inplace=True)
    with sys.stdout as fhout:
        resolved.to_csv(fhout, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--countsfile", type=str, required=True, help="Counts file of ASVs")
    parser.add_argument(
        "--clustfile",
        type=str, required=True,
        help="Taxonomy file for ASVs. Should also include a column with cluster designation. "
    )
    parser.add_argument(
        "--ranks",
        nargs="+",
        default=[
            "Kingdom",
            "Phylum",
            "Class",
            "Order",
            "Family",
            "Genus",
            "Species",
            "BOLD_bin",
        ],
        help="Ranks to include in the output (default: Kingdom Phylum Class Order Family Genus Species BOLD_bin))",
    )
    parser.add_argument(
        "--clust_column",
        type=str,
        default="cluster",
        help="Name of cluster column (default: 'cluster')"
    )
    parser.add_argument(
        "--consensus_threshold",
        type=int,
        default=80,
        help="Threshold (in %%) at which to assign taxonomy to a cluster (default: 80))",
    )
    parser.add_argument(
        "--consensus_ranks",
        nargs="+",
        default=["Family", "Genus", "Species", "BOLD_bin"],
        help="Ranks to use for calculating consensus. Must be present in the clustfile (default: Family Genus Species BOLD_bin))",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=10000,
        help="If countsfile is very large, specify chunksize to read it in a number of lines at a time",
    )
    parser.add_argument("--nrows", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    main(args)
