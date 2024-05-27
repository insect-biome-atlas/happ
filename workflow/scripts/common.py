#!/usr/bin/env python

from Bio import SeqIO
import pandas as pd
import gzip
import os
import shutil
import logging


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
    _filtered_ids = []
    # Read fasta file and write a new zipped fasta file with filtered seqs
    with gzip.open(outfile, "wt") as fhout:
        for seqid in filtered_ids:
            try:
                record = record_dict[seqid]
                fhout.write(f">{seqid}\n{record.seq}\n")
                _filtered_ids.append(seqid)
            except KeyError:
                continue
    return _filtered_ids


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
    os.makedirs(sm.params.tmpdir, exist_ok=True)
    logging.info(f"Reading taxonomic info from {sm.input.tax[0]}")
    taxdf = pd.read_csv(sm.input.tax[0], sep="\t", index_col=0, header=0)
    tax = sm.wildcards.tax
    split_rank = sm.params.split_rank
    dataf = taxdf.loc[taxdf[split_rank] == tax]
    logging.info(f"Indexing {sm.input.fasta}")
    record_dict = SeqIO.index(sm.input.fasta, "fasta")
    if ";size=" in list(record_dict.keys())[0]:
        record_dict = {
            k.split(";")[0] if ";" in k else k: v for k, v in record_dict.items()
        }
    logging.info(f"Writing sequences to {sm.params.fasta}")
    filtered_ids = write_fasta(record_dict, sm.params.fasta, list(dataf.index))
    logging.info(
        f"Extracting counts from {sm.input.counts[0]} and writing to {sm.params.total_counts} and {sm.params.counts}"
    )
    write_counts(
        countsfile=sm.input.counts[0],
        totalcounts=sm.params.total_counts,
        countsout=sm.params.counts,
        ids=filtered_ids,
    )
    logging.info(f"Moving {sm.params.total_counts} to {sm.output.total_counts}")
    shutil.move(sm.params.total_counts, sm.output.total_counts)
    logging.info(f"Moving {sm.params.fasta} to {sm.output.fasta}")
    shutil.move(sm.params.fasta, sm.output.fasta)
    logging.info(f"Moving {sm.params.counts} to {sm.output.counts}")
    shutil.move(sm.params.counts, sm.output.counts)
    shutil.rmtree(sm.params.tmpdir)


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
        for record in SeqIO.tmuxparse(fasta, "fasta"):
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
