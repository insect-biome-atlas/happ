#!/usr/bin/env python

from Bio.SeqIO import parse
import gzip


def format_swarm(sm):
    """
    Reformats input fasta file for use with swarm

    :param sm:
    :return:
    """
    counts = {}
    # Read total counts from input file
    with open(sm.input.counts, 'r') as fhin:
        for line in fhin:
            key, value = line.rstrip().split("\t")
            try:
                counts[key] = int(value)
            except ValueError:
                continue
    # Open fasta file and append counts as its being read
    with gzip.open(sm.input.fasta, 'rt') as fhin, gzip.open(sm.output.fasta, 'wt') as fhout:
        for record in parse(fhin, "fasta"):
            try:
                new_rec = f"{record.id}_{counts[record.id]}"
            except KeyError:
                continue
            if counts[record.id] >0:
                fhout.write(f">{new_rec}\n{record.seq}\n")


def main(sm):
    toolbox = {'format_swarm': format_swarm}
    toolbox[sm.rule](sm)


if __name__ == '__main__':
    main(snakemake)