#!/usr/bin/env python

from Bio.SeqIO import parse
import gzip
import shutil
import os
import pandas as pd
import re


def format_swarm(sm):
    """
    Reformats input fasta file for use with swarm

    :param sm:
    :return:
    """
    os.makedirs(sm.params.tmpdir, exist_ok=True)
    counts = {}
    # Read total counts from input file
    with open(sm.input.counts, "r") as fhin:
        for line in fhin:
            key, value = line.rstrip().split("\t")
            try:
                counts[key] = int(value)
            except ValueError:
                continue
    # Store sequences from fasta
    seqs = {}
    with gzip.open(sm.input.fasta, "rt") as fhin, gzip.open(
        sm.params.fasta, "wt"
    ) as fhout:
        for record in parse(fhin, "fasta"):
            try:
                new_rec = f"{record.id}_{counts[record.id]}"
                if str(record.seq) in seqs.keys():
                    seqs[str(record.seq)].append(new_rec)
                    continue
                else:
                    seqs[str(record.seq)] = [new_rec]
            except KeyError:
                continue
            if counts[record.id] > 0:
                fhout.write(f">{new_rec}\n{record.seq}\n")
    # Write dereplicated results (if any)
    with open(sm.output.derep, "w") as fhout:
        for k, v in seqs.items():
            fhout.write(f"{' '.join(v)}\n")
    shutil.move(sm.params.fasta, sm.output.fasta)
    shutil.rmtree(sm.params.tmpdir)


def get_cluster_members(f):
    col1 = []
    col2 = []
    regex = re.compile(r"_\d+$")
    with open(f, "r") as fhin:
        for i, line in enumerate(fhin, start=1):
            items = line.rstrip().rsplit()
            col1 += [regex.sub("", x) for x in items]
            col2 += [f"cluster{i}"] * len(items)
    dataf = pd.DataFrame(data={"ASV": col1, "cluster": col2})
    return dataf.set_index("ASV")


def read_derep(f, dataf):
    derep = {}
    regex = re.compile(r"_\d+$")
    with open(f, "r") as fhin:
        for line in fhin:
            items = [regex.sub("", x) for x in line.rstrip().rsplit()]
            rep = items[0]
            clust = dataf.loc[rep, "cluster"]
            if len(items) > 1:
                for item in items[1:]:
                    derep[item] = clust
    derepdf = pd.DataFrame(derep, index=["cluster"]).T
    derepdf.index.name = "ASV"
    return derepdf


def swarm2tab(sm):
    tmpdir = os.path.expandvars(sm.params.tmpdir)
    tmpout = os.path.expandvars(sm.params.out)
    f = sm.input[0]
    f2 = sm.input[1]
    dataf = get_cluster_members(f)
    derepdf = read_derep(f2, dataf)
    dataf = pd.concat([dataf, derepdf])
    os.makedirs(tmpdir, exist_ok=True)
    dataf.to_csv(tmpout, sep="\t")
    shutil.move(tmpout, sm.output[0])
    os.removedirs(tmpdir)


def main(sm):
    toolbox = {"format_swarm": format_swarm, "swarm2tab": swarm2tab}
    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)
