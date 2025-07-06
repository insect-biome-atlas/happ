import pandas as pd
from argparse import ArgumentParser


def concat_files(files, has_header=True):
    """
    Concatenates multiple tab-separated value (TSV) files into a single DataFrame.

    Args:
        files (list of str): List of file paths to the TSV files to be concatenated.

    Returns:
        pandas.DataFrame: A DataFrame containing the concatenated data from all input files.

    The function reads each file into a DataFrame, ensuring that all DataFrames have the same columns
    by using the columns from the first file. It then concatenates all DataFrames along the row axis.
    """
    df = pd.DataFrame()
    if has_header:
        header = 0
    else:
        header = None
    for i, f in enumerate(files):
        _df = pd.read_csv(f, sep="\t", index_col=0, header=header)
        if has_header:
            if i == 0:
                cols = _df.columns
            else:
                _df = _df.loc[:, cols]
        df = pd.concat([df, _df], axis=0)
    return df


def main(args):
    # filter cluster taxonomy
    taxdf = pd.read_csv(args.taxonomy, sep="\t", index_col=0)
    clust_col = False
    if "cluster" in taxdf.columns:
        clust_col = True
        index_name = taxdf.index.name
        taxdf = taxdf.reset_index().set_index("cluster")
    retained = concat_files(args.retained)
    if retained.index.name == None:
        retained.index.name = "cluster"
    if args.singles:
        singles = pd.read_csv(
            args.singles, sep="\t", index_col=0, header=0, 
        )
    if len(singles) > 0:
        singles = singles.loc[:, retained.columns]
        retained = pd.concat([retained, singles])
    retained = taxdf.loc[retained.index]
    discarded_seqs = list(set(taxdf.index).difference(retained.index))
    discarded = taxdf.loc[discarded_seqs]
    # filter consensus taxonomy
    if args.consensus_taxonomy:
        cons_tax = pd.read_csv(args.consensus_taxonomy, sep="\t", index_col=0)
        cons_tax.loc[retained.index.unique()].to_csv(
            f"{args.outdir}/noise_filtered_cluster_consensus_taxonomy.tsv", sep="\t"
        )
    # filter counts
    if args.counts:
        counts = pd.read_csv(args.counts, sep="\t", index_col=0)
        counts.loc[retained.index.unique()].to_csv(
            f"{args.outdir}/noise_filtered_cluster_counts.tsv", sep="\t"
        )
    if clust_col:
        retained = retained.reset_index().set_index(index_name)
        discarded = discarded.reset_index().set_index(index_name)
    # write discarded
    discarded.to_csv(
        f"{args.outdir}/discarded_cluster_taxonomy.tsv", sep="\t"
    )
    # write retained    
    retained.to_csv(
        f"{args.outdir}/noise_filtered_cluster_taxonomy.tsv", sep="\t"
    )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-r",
        "--retained",
        nargs="+",
        help="One or more files with retained sequences",
        required=True,
    )
    parser.add_argument(
        "-t",
        "--taxonomy",
        type=str,
        help="TSV file with taxonomic and cluster/OTU assignments",
        required=True,
    )
    parser.add_argument(
        "-o", "--outdir", type=str, help="Output directory for files", required=True
    )
    parser.add_argument(
        "-c",
        "--counts",
        type=str,
        help="Counts file for OTU/clusters",
    )
    parser.add_argument(
        "--consensus_taxonomy",
        type=str,
        help="TSV file with consensus taxonomy for clusters/OTUs",
    )
    parser.add_argument(
        "-s",
        "--singles",
        type=str,
        help="TSV file with single clusters/OTUs. Sequence ids in first column and cluster/OTU name in second.",
    )
    args = parser.parse_args()
    main(args)
