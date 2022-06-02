#!/usr/bin/env python

from Bio.SeqIO import parse
import gzip


def sum_counts(f, ids=[]):
    """
    Sum total counts of ASVs across all samples
    :param f: Input ASV table
    :param ids: Limit summing to optional list of ids
    :return: dictionary with summed counts
    """
    counts = {}
    with open(f, 'r') as fhin:
        for i, line in enumerate(fhin):
            if i == 0:
                continue
            asv = line.rsplit()[0]
            if len(ids)>0:
                if asv not in ids:
                    continue
            summed_count = sum([int(x) for x in line.rsplit()[1:]])
            counts[asv] = summed_count
    return counts


def filter_seqs(sm):
    """
    Create a table with total counts per ASV and write a new fasta file
    containing ASVs with >0 total counts across samples.

    :param sm:
    :return:
    """
    counts = sum_counts(sm.input.counts[0])
    filtered_ids = []
    # Write total counts to file
    with open(sm.output.total_counts, 'w') as fhout:
        fhout.write("Representative_Sequence\ttotal\n")
        for seqid, count in counts.items():
            if count > 0:
                fhout.write(f"{seqid}\t{count}\n")
                filtered_ids.append(seqid)
    # Read fasta file and write a new zipped fasta file with filtered seqs
    with open(sm.input.fasta[0], 'r') as fhin, gzip.open(sm.output.fasta, 'wt') as fhout:
        for record in parse(fhin, "fasta"):
            if record.id in filtered_ids:
                fhout.write(f">{record.id}\n{record.seq}\n")


def main(sm):
    toolbox = {'filter_seqs': filter_seqs}
    toolbox[sm.rule](sm)


if __name__ == '__main__':
    main(snakemake)