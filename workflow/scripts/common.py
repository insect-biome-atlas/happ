#!/usr/bin/env python

from Bio.SeqIO import parse
import pandas as pd
import gzip
import os
import shutil


def read_taxa(config):
    rundir = config["rundir"]
    split_rank = config["split_rank"]
    # See if the list of taxa already exists
    if os.path.exists(f"data/{rundir}/{split_rank}.txt"):
        taxa = []
        with open(f"data/{rundir}/{split_rank}.txt", 'r') as fhin:
            for line in fhin:
                taxa.append(line.rstrip())
        return taxa
    # Set file paths
    counts_file = f"data/{rundir}/asv_counts.tsv"
    fasta_file = f"data/{rundir}/asv_seqs.fasta"
    tax_file = f"data/{rundir}/asv_taxa.tsv"
    # Read taxonomy dataframe
    tax_df = pd.read_csv(tax_file, sep="\t", index_col=0)
    # Read counts and sum
    total_counts = sum_counts(counts_file)
    dataf = pd.DataFrame(total_counts, index=["count"]).T
    # Filter dataframe to asvs with sum >0
    dataf = dataf.loc[dataf["count"]>0]
    filtered_ids = []
    # Go through fasta file and make sure to filter to ids present there
    with open(fasta_file, 'r') as fhin:
        for record in parse(fhin, "fasta"):
            if record.id in dataf.index:
                filtered_ids.append(record.id)
    dataf = tax_df.loc[filtered_ids]
    # Count size of remaining taxa
    rank_size = dataf.groupby(split_rank).size()
    # Filter to taxa with at least 2 sequences
    filtered_taxa = list(rank_size.loc[rank_size>1].index)
    dataf = dataf.loc[dataf[split_rank].isin(filtered_taxa)]
    taxa = list(dataf[split_rank].unique())
    if not os.path.exists(f"data/{rundir}/{split_rank}.txt"):
        with open(f"data/{rundir}/{split_rank}.txt", 'w') as fhout:
            for tax in taxa:
                fhout.write(f"{tax}\n")
    return taxa


def sum_counts(f, fhout=None, sum_counts=True, ids=None):
    """
    Sum total counts of ASVs across all samples
    :param f: Input ASV table
    :param fhout: optional output file handle
    :param sum_counts: whether to do any summing or not
    :param ids: Limit summing to optional list of ids
    :return: dictionary with summed counts
    """
    if ids is None:
        ids = []
    counts = {}
    with open(f, 'r') as fhin:
        for i, line in enumerate(fhin):
            if i == 0:
                if fhout is not None:
                    fhout.write(line)
                continue
            asv = line.rsplit()[0]
            if len(ids)>0:
                if asv not in ids:
                    continue
            if sum_counts:
                summed_count = sum([int(x) for x in line.rsplit()[1:]])
                counts[asv] = summed_count
            if fhout is not None:
                fhout.write(line)
    return counts


def write_total(total_counts, outfile, ids=None):
    """
    Write a two column table with ASV id in first column and total counts in second column
    :param total_counts: dictionary of summed up counts for ASVs
    :param outfile: output file
    :param ids: optional ids to limit writing
    :return: list of ids with total_count > 0
    """
    if ids is None:
        ids = []
    filtered_ids = []
    if len(ids) > 0:
        items = ids
    else:
        items = list(total_counts.keys())
    # Write total counts to file
    with open(outfile, 'w') as fhout:
        fhout.write("Representative_Sequence\ttotal\n")
        for seqid in items:
            try:
                count = total_counts[seqid]
            except KeyError:
                continue
            if count > 0:
                fhout.write(f"{seqid}\t{count}\n")
                filtered_ids.append(seqid)
    return filtered_ids


def write_fasta(infile, outfile, filtered_ids):
    """
    Read a fasta file and write sequences to outfile if present in filtered_ids list

    :param infile: input fasta file
    :param outfile: output fasta file (gzipped)
    :param filtered_ids: list of ids to write
    :return: new filtered list
    """
    _filtered_ids = []
    # Read fasta file and write a new zipped fasta file with filtered seqs
    with open(infile, 'r') as fhin, gzip.open(outfile,'wt') as fhout_fasta:
        for record in parse(fhin, "fasta"):
            if record.id in filtered_ids:
                fhout_fasta.write(f">{record.id}\n{record.seq}\n")
                _filtered_ids.append(record.id)
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
    os.makedirs(sm.params.tmpdir, exist_ok=True)
    total_counts = sum_counts(sm.input.counts[0])
    taxdf = pd.read_csv(sm.input.tax[0], sep="\t", index_col=0, header=0)
    tax = sm.wildcards.tax
    split_rank = sm.params.split_rank
    dataf = taxdf.loc[taxdf[split_rank] == tax]
    filtered_ids = write_total(total_counts, sm.params.total_counts, dataf.index)
    filtered_ids = write_fasta(sm.input.fasta[0], sm.params.fasta, filtered_ids)
    with gzip.open(sm.params.counts, 'wt') as fhout:
        _ = sum_counts(sm.input.counts[0], fhout=fhout, sum_counts=False, ids=filtered_ids)
    shutil.move(sm.params.total_counts, sm.output.total_counts)
    shutil.move(sm.params.fasta, sm.output.fasta)
    shutil.move(sm.params.counts, sm.output.counts)
    shutil.rmtree(sm.params.tmpdir)


def collate(sm):
    cluster_num = 0
    prog = sm.wildcards.prog
    with open(sm.output[0], 'w') as fhout:
        fhout.write("asv\tcluster\ttax\torg_cluster\n")
        for f in sm.input:
            cluster_map = {}
            tax = os.path.basename(os.path.dirname(f))
            cluster_map[tax] = {}
            with open(f, 'r') as fhin:
                for i, line in enumerate(fhin):
                    if i==0:
                        continue
                    asv, cluster = line.rstrip().rsplit()
                    try:
                        newcluster = cluster_map[tax][cluster]
                    except KeyError:
                        newcluster = f"cluster{cluster_num}"
                        cluster_map[tax][cluster] = newcluster
                        cluster_num+=1
                    fhout.write(f"{asv}\t{newcluster}\t{tax}\t{cluster}\n")



def main(sm):
    toolbox = {'filter_seqs': filter_seqs,
               'collate': collate}
    toolbox[sm.rule](sm)


if __name__ == '__main__':
    main(snakemake)