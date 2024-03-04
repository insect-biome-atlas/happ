 #!/usr/bin/env python

import pandas as pd
from argparse import ArgumentParser
import sys
from Bio.SeqIO import parse
from tqdm import tqdm

def read_fasta(fasta):
    """
    Read fasta file and return ids
    """
    with open(fasta, "r") as f:
        return [record.id for record in parse(fasta, "fasta")]

def generate_reader(f, chunksize):
    """
    Sets up a reader with pandas. Handles both chunksize>=1 and chunksize=None

    :param f: Input file
    :param chunksize: Number of rows to read per chunk
    :return:
    """
    r = pd.read_csv(f, sep="\t", index_col=0, header=0, chunksize=chunksize)
    if chunksize is not None:
        return r
    return [r]

def read_counts(countsfile, n_asvs):
    reader = generate_reader(countsfile, chunksize=10000)
    asv_sum = pd.DataFrame()
    asv_obs = pd.DataFrame()
    for i, df in enumerate(tqdm(reader,unit=" chunks", desc="Reading counts", total=round(n_asvs / 10000))):
        _asv_sum = pd.DataFrame(df.sum(axis=1), columns=["sum"])
        asv_sum = pd.concat([asv_sum, _asv_sum])
        _asv_obs = pd.DataFrame(df.gt(0).sum(axis=1), columns=["n_obs"])
        asv_obs = pd.concat([asv_obs, _asv_obs])
    return pd.merge(asv_sum, asv_obs, left_index=True, right_index=True)

def count_clusters(df):
    return len(df["cluster"].unique())

def count_species(df):
    return len(df["Species"].unique())

def eval_chimeras(df):
    out = {}
    out["total_reads"] = df["sum"].sum()
    out["total_asvs"] = df.shape[0]
    out["trusted_asvs"] = df.loc[df.trusted == True].shape[0]
    chimera_df = df.loc[df.chimera == True]
    nonchimera_df = df.loc[df.nonchimera == True]
    out["total_clusters"] = len(nonchimera_df["cluster"].unique())
    trusted_chimeras = chimera_df.loc[(chimera_df.trusted == True)]
    trusted_removed = trusted_chimeras.shape[0]
    out["trusted_ASVs_removed"] = trusted_removed
    trusted_reads_removed = trusted_chimeras["sum"].sum()
    out["trusted_reads_removed"] = trusted_reads_removed
    asvs_removed = chimera_df.shape[0]
    out["asvs_removed"] = asvs_removed
    reads_removed = chimera_df["sum"].sum()
    out["reads_removed"] = reads_removed
    out["single_sample_asv_removed"] = len(chimera_df.loc[chimera_df["n_obs"] == 1])
    out["single_sample_asvs_remaining"] = len(nonchimera_df.loc[nonchimera_df["n_obs"] == 1])
    classified_df = nonchimera_df.loc[~nonchimera_df["Species"].str.startswith("unclassified")]
    cluster_species_ratio = len(classified_df["cluster"].unique()) / len(classified_df["Species"].unique())
    out["cluster_species_ratio"] = cluster_species_ratio
    species_per_cluster = (
        classified_df.groupby("cluster").apply(count_species)
    )
    multi_species_clusters = species_per_cluster.loc[species_per_cluster > 1].shape[0]
    out["multi_species_clusters"] = multi_species_clusters
    clusters_per_species = (
        classified_df.groupby("Species").apply(count_clusters)
    )
    multi_cluster_species = clusters_per_species.loc[clusters_per_species > 1].shape[0]
    out["multi_cluster_species"] = multi_cluster_species
    return pd.DataFrame(out, index=["eval"]).T

def main(args):
    trusted = []
    sys.stderr.write("#Reading taxonomies\n")
    tax = pd.read_csv(args.taxfile, sep="\t", index_col=0)
    clust = pd.read_csv(args.cluster_taxonomy, sep="\t", index_col=0, usecols=[0,1])
    df = pd.merge(tax, clust, left_index=True, right_index=True, how="left")
    df = df.assign(trusted=pd.Series([False]*df.shape[0], index=df.index))
    df = df.assign(chimera=pd.Series([False]*df.shape[0], index=df.index))
    df = df.assign(nonchimera=pd.Series([False]*df.shape[0], index=df.index))
    if args.trusted_fasta:
        trusted = read_fasta(args.trusted_fasta)
        df.loc[trusted, "trusted"] = True
    chimeras = read_fasta(args.chimera_fasta)
    df.loc[chimeras, "chimera"] = True
    nonchimeras = read_fasta(args.nonchimera_fasta)
    df.loc[nonchimeras, "nonchimera"] = True
    counts = read_counts(args.counts, df.shape[0])
    df = pd.merge(df, counts, left_index=True, right_index=True)
    out = eval_chimeras(df)
    out.to_csv(sys.stdout, sep="\t")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-t", "--taxfile", help="ASV taxonomy file", required=True)
    parser.add_argument("-c", "--cluster_taxonomy", help="Cluster taxonomy file", required=True)
    parser.add_argument("--chimera_fasta", help="Chimera fasta file", required=True)
    parser.add_argument("--nonchimera_fasta", help="Non-chimera fasta file", required=True)
    parser.add_argument("--counts", help="Counts file of ASVs", required=True)
    parser.add_argument("--trusted_fasta", help="Trusted fasta file")
    args = parser.parse_args()
    main(args)