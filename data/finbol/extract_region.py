#!/usr/bin/env python
###################################
from Bio.SeqIO import parse

infile = "finbol.msa.nohmm.fasta"
outfile = "finbol.asv-region.fasta"
with open(outfile, 'w') as fhout:
    for record in parse(infile, 'fasta'):
        seq = record.seq[283:750]
        if seq.startswith("-"):
            continue
        fhout.write(f">{record.description}\n{seq}\n")

