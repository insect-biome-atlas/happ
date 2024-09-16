#!/usr/bin/env python


"""
    Output should look like this:
    seqid1	k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Legionellales; f__Legionellaceae; g__Legionella; s__
    seqid2	k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae; g__Flavobacterium; s__

    fasta output should be:
    >seqid1
    ATGCGGGCTAGAGTAGCGAT...
"""

import pandas as pd
from Bio.SeqIO import parse
from argparse import ArgumentParser


def read_seqs(fasta):
    seqs = {}
    tax = {}
    for record in parse(fasta, "fasta"):
        seqid, desc = (record.id).split(";")
        desc = desc.replace("tax=", "")
        items = desc.split(",")
        lineage = "; ".join([x[0]+"__"+x[2:] for x in items])
        tax[seqid] = lineage
        seqs[seqid] = record
    return seqs, tax

def main(args):
    fasta = args.fasta
    tsv_output = args.tsv_output
    fasta_output = args.fasta_output
    seqs, tax = read_seqs(fasta)
    data = pd.DataFrame(tax, index=["Taxon"]).T
    data.index.name = "Feature ID"
    data.to_csv(tsv_output, sep="\t")
    with open(fasta_output, 'w') as fhout:
        for seqid, record in seqs.items():
            fhout.write(f">{seqid}\n{str(record.seq)}\n")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("fasta", help="Fasta file with sequences")
    parser.add_argument("tsv_output", help="Output file name")
    parser.add_argument("fasta_output", help="Output file name")
    args = parser.parse_args()
    main(args)