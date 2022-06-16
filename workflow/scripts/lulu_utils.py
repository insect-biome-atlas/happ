#!/usr/bin/env python

import pandas as pd
import os
import shutil


def lulu2tab(sm):
    df = pd.read_csv(sm.input[0], sep="\t", index_col=0)
    os.makedirs(sm.params.tmpdir)
    d = {}
    dataf = df.loc[df.curated=="merged"].sort_values("rank")
    for row in dataf.iterrows():
        asv, parent, rank = row[0], row[1]["parent_id"], row[1]["rank"]
        d[asv] = {"parent": parent, "rank": rank}
        d[parent] = {"parent": parent, "rank": rank}
    dataf = pd.DataFrame(d).T
    dataf.index.name = "asv"
    dataf = dataf.sort_values("rank")
    dataf = dataf.drop("rank", axis=1).reset_index()
    dataf = dataf.rename(index=lambda x: f"cluster{x}").drop("parent", axis=1)
    dataf.index.name = "cluster"
    dataf.to_csv(sm.params.out, sep="\t")
    shutil.move(sm.params.out, sm.output[0])
    os.rmdir(sm.params.tmpdir)


def main(sm):
    toolbox = {'lulu2tab': lulu2tab}
    toolbox[sm.rule][sm]


if __name__ == '__main__':
    main(snakemake)