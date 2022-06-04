#!/usr/bin/env python

from Bio.SeqIO import parse
import gzip


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


def write_total(total_counts, outfile):
    """
    Write a two column table with ASV id in first column and total counts in second column
    :param total_counts: dictionary of summed up counts for ASVs
    :param outfile: output file
    :return: list of ids with total_count > 0
    """
    filtered_ids = []
    # Write total counts to file
    with open(outfile, 'w') as fhout:
        fhout.write("Representative_Sequence\ttotal\n")
        for seqid, count in total_counts.items():
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
    total_counts = sum_counts(sm.input.counts[0])
    filtered_ids = write_total(total_counts, sm.output.total_counts)
    filtered_ids = write_fasta(sm.input.fasta[0], sm.output.fasta, filtered_ids)
    with gzip.open(sm.output.counts, 'wt') as fhout:
        _ = sum_counts(sm.input.counts[0], fhout=fhout, sum_counts=False, ids=filtered_ids)


def main(sm):
    toolbox = {'filter_seqs': filter_seqs}
    toolbox[sm.rule](sm)


if __name__ == '__main__':
    main(snakemake)