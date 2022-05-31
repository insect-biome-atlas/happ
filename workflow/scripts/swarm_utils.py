#!/usr/bin/env python

from Bio.SeqIO import parse
from common import sum_counts


def format_swarm(sm):
    counts = sum_counts(sm.input.counts[0])
    with open(sm.input.fasta[0], 'r') as fhin, open(sm.output.fasta, 'w') as fhout:
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