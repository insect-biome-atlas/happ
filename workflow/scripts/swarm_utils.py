#!/usr/bin/env python

from Bio.SeqIO import parse
import gzip
import shutil
import os
import pandas as pd


def format_swarm(sm):
    """
    Reformats input fasta file for use with swarm

    :param sm:
    :return:
    """
    os.makedirs(sm.params.tmpdir, exist_ok=True)
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
    with gzip.open(sm.input.fasta, 'rt') as fhin, gzip.open(sm.params.fasta, 'wt') as fhout:
        for record in parse(fhin, "fasta"):
            try:
                new_rec = f"{record.id}_{counts[record.id]}"
            except KeyError:
                continue
            if counts[record.id] >0:
                fhout.write(f">{new_rec}\n{record.seq}\n")
    shutil.move(sm.params.fasta, sm.output.fasta)
    shutil.rmtree(sm.params.tmpdir)


def get_cluster_members(f):
    col1 = []
    col2 = []
    with open(f, 'r') as fhin:
        for i, line in enumerate(fhin, start=1):
            items = line.rstrip().rsplit()
            col1 += [x.split("_")[0] for x in items]
            col2 += [f"cluster{i}"] * len(items)
    dataf = pd.DataFrame(data={'ASV': col1, 'cluster': col2})
    return dataf.set_index("ASV")


def swarm2tab(sm):
    f = sm.input[0]
    dataf = get_cluster_members(f)
    os.makedirs(sm.params.tmpdir, exist_ok=True)
    dataf.to_csv(sm.params.out, sep="\t")
    shutil.move(sm.params.out, sm.output[0])
    os.removedirs(sm.params.tmpdir)

def main(sm):
    toolbox = {'format_swarm': format_swarm,
               'swarm2tab': swarm2tab}
    toolbox[sm.rule](sm)


if __name__ == '__main__':
    main(snakemake)
