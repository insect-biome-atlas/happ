#!/usr/env/bin python

import argparse
import pandas as pd
import sys
import logging
import tqdm


def generate_reader(f, chunksize, nrows):
    """
    Sets up a reader with pandas. Handles both chunksize>=1 and chunksize=None

    :param f: Input file
    :param chunksize: Number of rows to read per chunk
    :param nrows: Number of total rows to read
    :return:
    """
    if nrows == 0:
        nrows = None
    r = pd.read_csv(
        f, sep="\t", index_col=0, header=0, nrows=nrows, chunksize=chunksize
    )
    if chunksize is not None:
        return r
    return [r]


def read_metadata(f, index_name="sampleID_NGI"):
    logging.info(f"Reading metadata from {f}")
    #sys.stderr.write("####\n" f"Reading metadata from {f}\n")
    df = pd.read_csv(f, sep="\t", header=0, comment="#")
    return df.set_index(index_name)


def read_counts(
    countsfile,
    metadata,
    split_col,
    blanks,
    chunksize,
    nrows=None,
):
    """
    Read the counts file in chunks, if list of blanks is given, count occurrence
    in blanks and return as a column <in_n_blanks>. Also calculate max and
    sum for each ASV.
    """
    reader = generate_reader(countsfile, chunksize=chunksize, nrows=nrows)
    logging.info(f"Reading counts from {countsfile}")
    #sys.stderr.write("####\n" f"Reading counts from {countsfile}\n")
    data = {}
    warnings = []
    n_asvs = n_samples = 0
    n_datasets = []
    for i, df in enumerate(
        tqdm.tqdm(
            reader,
            unit=" chunks",
            leave=False,
        )
    ):
        n_asvs += df.shape[0]
        if i == 0:
            n_samples = df.shape[1]
            sample_names = list(df.columns)
            if not split_col:
                split_col = "dataset"
                metadata[split_col] = [split_col]*metadata.shape[0]
        # get unique values of split_col
        split_col_vals = metadata[split_col].unique()
        for val in split_col_vals:
            if val not in data.keys():
                data[val] = pd.DataFrame()
            # get samples corresponding to split
            val_samples = metadata.loc[metadata[split_col] == val].index
            # get intersection of val_samples and the df columns
            val_samples_intersect = list(set(val_samples).intersection(df.columns))
            if len(val_samples_intersect) == 0 and i == 0:
                warnings.append(
                    "####\n" f"No samples found in counts data for {val}, skipping...\n"
                )
                continue
            n_datasets.append(val)
            # get samples in val_samples missing from df columns
            missing_samples = set(val_samples).difference(val_samples_intersect)
            if len(missing_samples) > 0 and i == 0:
                warnings.append(
                    "####\n"
                    f"WARNING: {len(missing_samples)} samples in metadata file are missing from counts file for {val}:\n"
                )
                warnings.append(f"{', '.join(missing_samples)} \n")
            # split the counts dataframe
            split_df = df.loc[:, val_samples_intersect]
            # get blanks for this dataset
            val_blanks = list(set(blanks).intersection(val_samples_intersect))
            # get blanks in metadata actually present in counts data
            val_blanks_intersect = list(
                set(val_blanks).intersection(val_samples_intersect)
            )
            # get blanks missing
            # missing_blanks = set(val_blanks).difference(val_blanks_intersect)
            if i == 0:
                warnings.append(
                    "####\n" f"{len(val_blanks_intersect)} blanks found for {val}\n"
                )
            # calculate ASV sum (remove blanks)
            asv_sum = pd.DataFrame(
                split_df.drop(val_blanks_intersect, axis=1).sum(axis=1),
                columns=["ASV_sum"],
            )
            # calculate ASV max
            asv_max = pd.DataFrame(
                split_df.drop(val_blanks_intersect, axis=1).max(axis=1),
                columns=["ASV_max"],
            )
            _dataframe = pd.merge(asv_sum, asv_max, left_index=True, right_index=True)
            if len(val_blanks_intersect) > 0:
                blank_counts = split_df.loc[:, val_blanks_intersect]
                # calculate occurrence in blanks
                asv_blank_count = pd.DataFrame(
                    blank_counts.gt(0).sum(axis=1), columns=["in_n_blanks"]
                )
                asv_blank_count["in_percent_blanks"] = (
                    asv_blank_count.div(len(val_blanks_intersect)) * 100
                )
                _dataframe = pd.merge(
                    _dataframe, asv_blank_count, left_index=True, right_index=True
                )
            data[val] = pd.concat([data[val], _dataframe])
    logging.info(f"Read counts for {n_asvs} ASVs in {n_samples} samples and {len(set(n_datasets))} datasets")
    return data

def main(args):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(message)s",
    )
    # Read metadata
    metadata = read_metadata(args.metadata, index_name=args.sample_id_col)
    logging.info(f"Read metadata with {metadata.shape[0]} samples")
    # Extract blanks from metadata
    blanks = metadata.loc[metadata[args.sample_type_col].isin(args.blank_val)].index.tolist()
    logging.info(f"Found {len(blanks)} blanks in metadata")
    # Read cluster assignments
    logging.info(f"Reading cluster assignments from {args.clustfile}")
    clustdf = pd.read_csv(args.clustfile, sep="\t", index_col=0)
    # Read counts (returns a dictionary where keys are datasets (or equivalent) and values are dataframes)
    # the dataframes have the columns ASV_sum, ASV_max, in_n_blanks, in_percent_blanks, for example:
    #                                   ASV_sum  ASV_max  in_n_blanks  in_percent_blanks
    # ASV_ID                                                                            
    # c4a912c251d80369193639cd88e5a82f      121       67            1              100.0
    counts = read_counts(
        countsfile=args.countsfile,
        metadata=metadata,
        split_col=args.split_col,
        blanks=blanks,
        chunksize=args.chunksize,
        nrows=args.nrows,
    )
    # Iterate over datasets and store ASVs to remove
    blank_asvs = []
    for dataset, dataframe in counts.items():
        blank_asvs += dataframe.loc[dataframe.in_percent_blanks>args.max_blank_occurrence].index.tolist()
    blank_asvs = list(set(blank_asvs))
    logging.info(f"Found {len(blank_asvs)} ASVs present in >{args.max_blank_occurrence}% of blanks")
    # Get clusters to remove
    blank_asvs_intersect = list(set(blank_asvs).intersection(clustdf.index))
    logging.info(f"Found {len(blank_asvs_intersect)}/{len(blank_asvs)} ASVs in clusters to remove")
    blank_clusters = clustdf.loc[blank_asvs_intersect, args.clustcol].unique().tolist()
    # Get ASVs in blank_clusters
    asvs_in_blank_clusters = clustdf.loc[clustdf[args.clustcol].isin(blank_clusters)].index.tolist()
    logging.info(f"Will remove {len(blank_clusters)} clusters with {len(asvs_in_blank_clusters)} ASVs")
    if args.outfile:
        logging.info(f"Writing filtered output to {args.outfile}")
        clustdf.drop(asvs_in_blank_clusters).to_csv(args.outfile, sep="\t")
    else:
        logging.info("Writing filtered output to stdout")
        clustdf.drop(asvs_in_blank_clusters).to_csv(sys.stdout, sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        """
        This script removes ASVs that belong to clusters where at least
        one ASV occurs in > max_blank_occurrence% of blanks.
        """
    )
    io_group = parser.add_argument_group("input/output")
    io_group.add_argument("--clustfile", type=str, help="Tab-separated file with cluster designations for ASVs", required=True)
    io_group.add_argument("--countsfile", type=str, help="Counts file of ASVs", required=True)
    io_group.add_argument("--outfile", type=str, help="Filtered counts file")
    io_group.add_argument("--metadata", type=str, help="Metadata file with sample information", required=True)
    io_group.add_argument(
        "--sample_id_col",
        type=str,
        help="Name of column in metadata file that contains sample ids",
    )
    io_group.add_argument(
        "--split_col",
        type=str,
        help="Name of column in metadata file by which to split samples by prior to cleaning by blanks",
    )
    io_group.add_argument(
        "--clustcol",
        type=str,
        help="Column in clustfile with cluster assignments (default: 'cluster')",
        default="cluster",
    )
    params_group = parser.add_argument_group("params")
    params_group.add_argument(
        "--sample_type_col",
        type=str,
        default="lab_sample_type",
        help="Use this column in metadata to identify sample type (default: 'lab_sample_type')",
    )
    params_group.add_argument(
        "--blank_val",
        type=str,
        nargs="+",
        default=["buffer_blank", "extraction_neg", "pcr_neg"],
        help="Values in <sample_type_col> that identify blanks (default: ['buffer_blank', 'extraction_neg', 'pcr_neg'])",
    )
    params_group.add_argument(
        "--max_blank_occurrence",
        type=int,
        help="Remove ASVs occurring in <max_blank_occurrence>%% of blank samples (default: 5)",
        default=5,
    )
    debug_group = parser.add_argument_group("debug")
    debug_group.add_argument(
        "--chunksize",
        type=int,
        help="Size of chunks (in lines) to read from countsfile (default: 10000)",
        default=10000
    )
    debug_group.add_argument(
        "--nrows",
        type=int,
        help=argparse.SUPPRESS,
    )
    args = parser.parse_args()
    main(args)