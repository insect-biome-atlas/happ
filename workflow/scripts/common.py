#!/usr/bin/env python

from Bio import SeqIO
import pandas as pd
import gzip
import os
import shutil
import logging
from tqdm import tqdm


def mem_allowed(wildcards, threads):
    return max(threads * 6400, 6400)


def read_taxa(config):
    rundir = config["rundir"]
    split_rank = config["split_rank"]
    chimera_run=config["chimera"]["run_name"]
    method=config["chimera"]["method"]
    algo=config["chimera"]["algorithm"]
    taxa = []
    # See if the list of taxa already exists
    # If there exists a taxfile from chimera filtering, only containing taxa with ASVs remaining
    # after filtering, use that
    # Example: If split_rank is set to 'Family', look for a file 'Family.txt' in the chimera filtered output
    # and use that if present
    if os.path.exists(f"results/chimera/{rundir}/filtered/{chimera_run}/{method}.{algo}/{split_rank}.txt"):
        with open(f"results/chimera/{rundir}/filtered/{chimera_run}/{method}.{algo}/{split_rank}.txt", "r") as fhin:
            for line in fhin:
                taxa.append(line.rstrip())
        return taxa
    elif os.path.exists(f"data/{rundir}/{split_rank}.txt"):
        with open(f"data/{rundir}/{split_rank}.txt", "r") as fhin:
            for line in fhin:
                taxa.append(line.rstrip())
        return taxa
    # Set file paths
    fasta_file = f"data/{rundir}/asv_seqs.fasta"
    tax_file = f"data/{rundir}/asv_taxa.tsv"
    # Read taxonomy dataframe
    tax_df = pd.read_csv(tax_file, sep="\t", index_col=0)
    filtered_ids = []
    # Go through fasta file and make sure to filter to ids present there
    with open(fasta_file, "r") as fhin:
        for record in SeqIO.parse(fhin, "fasta"):
            if record.id in tax_df.index:
                filtered_ids.append(record.id)
    dataf = tax_df.loc[filtered_ids]
    taxa = list(dataf[split_rank].unique())
    if not os.path.exists(f"data/{rundir}/{split_rank}.txt"):
        with open(f"data/{rundir}/{split_rank}.txt", "w") as fhout:
            for tax in taxa:
                fhout.write(f"{tax}\n")
    return taxa


def write_counts(countsfile, totalcounts, countsout, ids):
    n = len(ids)
    written = 0
    with open(countsfile, "r") as fhin, open(totalcounts, "w") as fh_total, gzip.open(
        countsout, "wt"
    ) as fh_counts:
        for i, line in enumerate(fhin):
            if i == 0:
                fh_total.write("Representative_Sequence\ttotal\n")
                fh_counts.write(line)
                continue
            items = line.rstrip().rsplit()
            if not items[0] in ids:
                continue
            s = int(sum([float(x) for x in items[1:]]))
            fh_total.write(f"{items[0]}\t{s}\n")
            fh_counts.write(line)
            written += 1
            if written == n:
                break


def write_fasta(record_dict, outfile, filtered_ids):
    """
    Read a fasta file and write sequences to outfile if present in filtered_ids list

    :param infile: input fasta file
    :param outfile: output fasta file (gzipped)
    :param filtered_ids: list of ids to write
    :return: new filtered list
    """
    records = []
    # Read fasta file and write a new zipped fasta file with filtered seqs
    for seqid in filtered_ids:
        try:
            records.append(record_dict[seqid])
        except KeyError:
            continue
    if len(records) > 0:
        with gzip.open(outfile, "wt") as fhout:
            SeqIO.write(records, fhout, "fasta")
    return [x.id for x in records]


def filter_seqs(sm):
    """
    This rule reads the raw counts table and fasta file and:
    1. sums up ASV counts across all samples and outputs this to a 'total_counts' file
       for use with opticlust and swarm
    2. writes new counts table and fasta file containing only sequences with
       total abundance > 0 and matching ASV ids
    :param sm:
    :return:
    """
    logging.basicConfig(
        filename=sm.log[0],
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s - %(message)s",
    )
    logging.info(f"Reading taxonomic info from {sm.input.tax[0]}")
    taxdf = pd.read_csv(sm.input.tax, sep="\t", index_col=0, header=0)
    outdir=sm.output[0]
    os.makedirs(outdir, exist_ok=True)
    split_rank = sm.params.split_rank
    taxa = taxdf[split_rank].unique()
    logging.info(f"Found {len(taxa)} unique taxa")
    logging.info(f"Indexing {sm.input.fasta}")
    record_dict = SeqIO.index(sm.input.fasta, "fasta")
    if ";size=" in list(record_dict.keys())[0]:
        record_dict = {
            k.split(";")[0] if ";" in k else k: v for k, v in record_dict.items()
        }
    taxa_w_seqs = []
    for tax in tqdm(taxa, desc="Extracting sequences for taxa"):
        os.makedirs(f"{outdir}/{tax}", exist_ok=True)
        dataf = taxdf.loc[taxdf[split_rank] == tax]
        fasta = f"{outdir}/{tax}/asv_seqs.fasta.gz"
        logging.info(f"Writing sequences to {fasta}")
        filtered_ids = write_fasta(record_dict, fasta, list(dataf.index))
        if len(filtered_ids) == 0:
            continue
        taxa_w_seqs.append(tax)
        countsout = f"{outdir}/{tax}/asv_counts.tsv.gz"
        totalcounts = f"{outdir}/{tax}/total_counts.tsv"
        logging.info(
            f"Extracting counts from {sm.input.counts[0]} and writing to {totalcounts} and {countsout}"
        )
        write_counts(
            countsfile=sm.input.counts,
            totalcounts=totalcounts,
            countsout=countsout,
            ids=filtered_ids,
        )

def collate(sm):
    cluster_num = 0
    prog = sm.wildcards.prog
    with open(sm.output[0], "w") as fhout:
        fhout.write("asv\tcluster\ttax\torg_cluster\n")
        for f in sm.input:
            cluster_map = {}
            tax = os.path.basename(os.path.dirname(f))
            cluster_map[tax] = {}
            with open(f, "r") as fhin:
                for i, line in enumerate(fhin):
                    if i == 0:
                        continue
                    asv, cluster = line.rstrip().rsplit()
                    try:
                        newcluster = cluster_map[tax][cluster]
                    except KeyError:
                        newcluster = f"cluster{cluster_num}"
                        cluster_map[tax][cluster] = newcluster
                        cluster_num += 1
                    fhout.write(f"{asv}\t{newcluster}\t{tax}\t{cluster}\n")


def clean_fasta(sm):
    """
    Reads a tab-separated file with ASV ids in first column and a fasta file
    and outputs only the fasta sequences for ASVs in the tab-separated file.
    """
    fasta = sm.input.fasta
    tsv = sm.input.tsv
    logfile = sm.log[0]
    logging.basicConfig(
        level=logging.INFO,
        filename=logfile,
        filemode="w",
        format="%(asctime)s - %(message)s",
    )
    outfile = sm.output[0]
    df = pd.read_csv(tsv, sep="\t", index_col=0)
    with open(outfile, 'w') as fhout:
        for record in SeqIO.parse(fasta, "fasta"):
            if record.id in df.index:
                fhout.write(f">{record.description}\n{record.seq}\n")

def clean_counts(sm):
    """
    Reads a tab-separated file with ASV ids in first column and a counts file
    with samples as columns and ASV clusters as rows, and outputs only clusters found
    in the TSV file.
    """
    tsv = sm.input.tsv
    counts = sm.input.counts
    logfile = sm.log[0]
    logging.basicConfig(
        level=logging.INFO,
        filename=logfile,
        filemode="w",
        format="%(asctime)s - %(message)s",
    )
    df = pd.read_csv(tsv, sep="\t", index_col=0)
    clusters = df["cluster"].unique()
    with open(sm.input.counts, 'r') as fhin, open(sm.output[0], 'w') as fhout:
        for i, line in enumerate(fhin):
            if i==0:
                fhout.write(line)
                continue
            if line.rsplit()[0] in clusters:
                fhout.write(line)

def main(sm):
    toolbox = {"filter_seqs": filter_seqs, "collate": collate, "clean_counts": clean_counts, "clean_fasta": clean_fasta}
    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)
