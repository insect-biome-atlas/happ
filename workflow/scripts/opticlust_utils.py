#!/usr/bin/env python

from common import sum_counts
from Bio.SeqIO import parse


def mothur_sum(sm):
    counts = sum_counts(sm.input.counts[0])
    filtered_ids = []
    with open(sm.output.counts, 'w') as fhout:
        fhout.write("Representative_Sequence\ttotal\n")
        for seqid, count in counts.items():
            if count > 0:
                fhout.write(f"{seqid}\t{count}\n")
                filtered_ids.append(seqid)
    with open(sm.input.fasta[0], 'r') as fhin, open(sm.output.fasta, 'w') as fhout:
        for record in parse(fhin, "fasta"):
            if record.id in filtered_ids:
                fhout.write(f">{record.id}\n{record.seq}\n")


def main(sm):
    toolbox = {'mothur_sum': mothur_sum}
    toolbox[sm.rule](sm)


if __name__ == '__main__':
    main(snakemake)