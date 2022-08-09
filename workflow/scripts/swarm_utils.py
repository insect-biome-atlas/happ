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


def get_cluster_members(df):
    dataf = pd.DataFrame()
    for col in ["asv1", "asv2"]:
        _ = df.groupby([col, "cluster"]).first().reset_index().loc[:,
            [col, "cluster"]]
        _.columns = ["asv", "cluster"]
        dataf = pd.concat([dataf, _])
    dataf["cluster"] = [f"cluster{x}" for x in dataf["cluster"]]
    dataf.set_index("asv", inplace=True)
    return dataf


def swarm2tab(sm):
    df = pd.read_csv(sm.input[0], sep="\t", header=None,
                     names=["asv1", "asv2", "diffs", "cluster",
                            "cumulative_steps"])
    dataf = get_cluster_members(df)
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
