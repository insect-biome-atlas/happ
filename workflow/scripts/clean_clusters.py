#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd
import sys
import tqdm


def read_taxonomy(f):
    sys.stderr.write("####\n" f"Reading clusters and taxonomy from {f}\n")
    asv_taxa = pd.read_csv(f, sep="\t", index_col=0)
    sys.stderr.write(f"{asv_taxa.shape[0]} ASVs read\n")
    return asv_taxa


def read_counts(f, blanks, chunksize=10000, nrows=0):
    """
    Read the counts file in chunks, if list of blanks is given, count occurrence
    in blanks and return as a column <in_n_blanks>. Also calculate max and
    sum for each ASV.
    """
    if nrows == 0:
        nrows = None
    sys.stderr.write("####\n" f"Reading counts from {f}\n")
    dataframe = pd.DataFrame()
    for i, df in enumerate(
        tqdm.tqdm(
            pd.read_csv(f, sep="\t", index_col=0, chunksize=chunksize, nrows=nrows),
            unit=" chunks",
        )
    ):
        if i == 0:
            samples = df.shape[1]
        if len(blanks) > 0:
            blank_counts = df.loc[:, blanks]
            # calculate occurrence in blanks
            asv_blank_count = pd.DataFrame(
                blank_counts.gt(0).sum(axis=1), columns=["in_n_blanks"]
            )
        else:
            asv_blank_count = pd.DataFrame(
                data={"in_n_blanks": [0] * df.shape[0]}, index=df.index
            )
        # calculate ASV sum (remove blanks)
        asv_sum = pd.DataFrame(df.drop(blanks, axis=1).sum(axis=1), columns=["ASV_sum"])
        # calculate ASV max
        asv_max = pd.DataFrame(df.drop(blanks, axis=1).max(axis=1), columns=["ASV_max"])
        _dataframe = pd.merge(asv_sum, asv_max, left_index=True, right_index=True)
        _dataframe = pd.merge(
            _dataframe, asv_blank_count, left_index=True, right_index=True
        )
        dataframe = pd.concat([dataframe, _dataframe])
    sys.stderr.write(
        f"Read counts for {dataframe.shape[0]} ASVs in " f"{samples} samples\n"
    )
    return dataframe


def clean_by_taxonomy(dataframe, rank="Family"):
    """
    Removes ASVs if they are 'unassigned' at <rank> or <rank> contains '_X'
    """
    df = dataframe.copy()
    before = df.shape[0]
    sys.stderr.write("####\n" f"Removing ASVs unclassified at {rank}\n")
    cleaned = df.loc[
        (~df[rank].str.contains("_X+$")) & (~df[rank].str.startswith("unclassified"))
    ]
    after = cleaned.shape[0]
    sys.stderr.write(
        f"{before - after} ASVs removed, {cleaned.shape[0]} ASVs remaining\n"
    )
    return cleaned


def read_blanks(f):
    sys.stderr.write("####\n" f"Reading list of blanks from {f}\n")
    with open(f, "r") as fhin:
        blanks = [x.rstrip() for x in fhin.readlines()]
    sys.stderr.write(f"{len(blanks)} blanks read\n")
    return blanks


def clean_by_reads(dataframe, min_clust_count=3):
    """
    Remove clusters with a sum less than <min_reads> across samples
    """
    sys.stderr.write(
        "####\n" f"Removing ASVs in clusters with <{min_clust_count} total reads\n"
    )
    df = dataframe.copy()
    # Groupby cluster column and sum values in ASV_sum
    cl_sum = df.groupby("cluster").sum(numeric_only=True)
    # Get list of clusters to remove
    cl_to_remove = cl_sum.loc[cl_sum["ASV_sum"] < min_clust_count].index
    # Get list of ASVs in said clusters
    asvs_to_remove = df.loc[df["cluster"].isin(cl_to_remove)].index
    before = df.shape[0]
    df.drop(asvs_to_remove, inplace=True)
    after = df.shape[0]
    sys.stderr.write(f"{before - after} ASVs removed, {df.shape[0]} ASVs remaining\n")
    return df


def clean_by_blanks(dataframe, blanks, mode="asv", max_blank_occurrence=20):
    """
    Removes clusters with ASVs present in > <max_blank_occurrence>% of blanks
    """
    sys.stderr.write(
        "####\n" f"Removing {mode}s in >{max_blank_occurrence}% of blanks\n"
    )
    df = dataframe.copy()
    before = df.shape[0]
    # Calculate % occurrence in blanks
    df["in_percent_blanks"] = df["in_n_blanks"].div(len(blanks)) * 100
    # Get list of ASVs to remove
    to_remove = df.loc[df.in_percent_blanks > max_blank_occurrence].index
    if mode == "cluster":
        # Get clusters for said ASVs
        to_remove_cl = df.loc[to_remove, "cluster"]
        # Get all ASVs in said clusters
        to_remove = list(df.loc[df["cluster"].isin(list(to_remove_cl))].index)
    df.drop(to_remove, inplace=True)
    after = df.shape[0]
    sys.stderr.write(f"{before - after} ASVs removed, {df.shape[0]} ASVs remaining\n")
    return df


def main(args):
    # Read taxonomy + clusters
    asv_taxa = read_taxonomy(args.taxonomy)
    # Read blanks
    blanks = read_blanks(args.blanks) if args.blanks else []
    # Read counts
    counts = read_counts(
        f=args.counts, blanks=blanks, chunksize=args.chunksize, nrows=args.nrows
    )
    # Clean by taxonomy
    asv_taxa_cleaned = clean_by_taxonomy(dataframe=asv_taxa, rank=args.clean_rank)
    # Merge counts + taxonomy
    asv_taxa_cleaned = pd.merge(
        asv_taxa_cleaned, counts, left_index=True, right_index=True
    )
    # Clean by blanks
    if args.blanks:
        asv_taxa_cleaned = clean_by_blanks(
            dataframe=asv_taxa_cleaned,
            blanks=args.blanks,
            mode=args.blank_removal_mode,
            max_blank_occurrence=args.max_blank_occurrence,
        )
    # Clean by read sum
    asv_taxa_cleaned = clean_by_reads(
        dataframe=asv_taxa_cleaned, min_clust_count=args.min_clust_count
    )
    # Write to output
    with sys.stdout as fhout:
        sys.stderr.write(
            "####\n" f"Writing {asv_taxa_cleaned.shape[0]} ASVs to stdout\n"
        )
        asv_taxa.index.name = "ASV"
        asv_taxa_cleaned.to_csv(fhout, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser(
        """
        This script cleans clustering results by removing ASVs if:
        - unassigned or ambiguous taxonomic assignments (e.g.  
        'unclassified' or '_X' in rank labels) 
        - if belonging to clusters present in > max_blank_occurrence% of blanks
        - if belonging to clusters with < min_clust_count total reads
        """
    )
    io_group = parser.add_argument_group("input/output")
    io_group.add_argument("--counts", type=str, help="Counts file of ASVs")
    io_group.add_argument(
        "--taxonomy",
        type=str,
        help="Taxonomy file for ASVs. Should also "
        "include a"
        "column with cluster designation.",
    )
    io_group.add_argument(
        "--blanks", type=str, help="File with samples that are 'blanks'"
    )
    io_group.add_argument("--output", type=str, help="Output file with cleaned results")
    params_group = parser.add_argument_group("params")
    params_group.add_argument(
        "--clean_rank",
        type=str,
        default="Family",
        help="Remove ASVs unassigned at this taxonomic " "rank (default Family)",
    )
    params_group.add_argument(
        "--max_blank_occurrence",
        type=float,
        default=20,
        help="Remove ASVs occurring in clusters where at "
        "least one member is present in "
        "<max_blank_occurrence>%% of blank samples. "
        "(default 20)",
    )
    params_group.add_argument(
        "--blank_removal_mode",
        type=str,
        choices=["cluster", "asv"],
        default="asv",
        help="How to remove sequences based on "
        "occurrence in blanks. If 'asv' ("
        "default) remove "
        "only ASVs that occur in more than "
        "<max_blank_occurrence>%% of blanks. If "
        "'cluster', remove ASVs in clusters where "
        "one or more ASVs is above the "
        "<max_blank_occurrence> threshold",
    )
    params_group.add_argument(
        "--min_clust_count",
        type=int,
        default=3,
        help="Remove clusters with < <min_clust_count> "
        "summed across samples (default 3)",
    )
    debug_group = parser.add_argument_group("debug")
    debug_group.add_argument(
        "--chunksize",
        type=int,
        default=10000,
        help="Size of chunks (in lines) to read from " "countsfile",
    )
    debug_group.add_argument(
        "--nrows",
        type=int,
        default=0,
        help="Number of rows to read from countsfile (" "for testing purposes only)",
    )

    args = parser.parse_args()
    main(args)
