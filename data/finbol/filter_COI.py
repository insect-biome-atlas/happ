#!/usr/bin/env python
###################################
from Bio.SeqIO import parse

with open("finbol.COI.fasta", 'w') as fhout:
    for record in parse("fasta.fas", 'fasta'):
        desc = record.description
        gene = desc.split("|")[2]
        if gene == "COI-5P":
            fhout.write(f">{desc}\n{record.seq}\n")

