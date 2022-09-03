#!/usr/bin/env python

import pandas as pd
import os
import shutil


def get_cluster_members(df):
    """
    Create a dataframe of ASV and cluster ids
    Input is a table as such:

            total   spread 	parent_id 	curated 	rank
    asv1 	4014184 424 	asv1 	    parent 	    1
    asv2 	166839 	154 	asv2 	    parent 	    2
    asv3 	28 	    1 	    asv2 	    merged 	    155
    asv4    858     1       asv1        merged      1
    asv5    99      33      asv5        parent      100

    This is parsed into:
            parent 	cluster
    asv4 	asv1 	cluster0
    asv1 	asv1 	cluster0
    asv3 	asv2 	cluster1
    asv2 	asv2 	cluster1
    """
    d = {}
    dataf = df.loc[df.curated == "merged"].sort_values("rank")
    # Loop through dataframe and store parent rank
    for row in dataf.iterrows():
        asv, parent, rank = row[0], row[1]["parent_id"], row[1]["rank"]
        d[asv] = {"parent": parent}
        d[parent] = {"parent": parent}
    # Create a new dataframe with asv id -> parent id
    dataf = pd.DataFrame(d).T
    dataf.index.name = "asv"
    # Create a dictionary with unique parent ids as keys and cluster ids as values
    cluster_nums = dict(zip(dataf.parent.unique(), [f"cluster{i}" for i in range(0,len(dataf.parent.unique())+1)]))
    # Make dictionary into a dataframe and merge with asv ids
    cluster_df = pd.DataFrame(cluster_nums, index=["cluster"]).T
    dataf = pd.merge(dataf, cluster_df, left_on="parent", right_index=True)
    return dataf


def lulu2tab(sm):
    df = pd.read_csv(sm.input[0], sep="\t", index_col=0)
    os.makedirs(sm.params.tmpdir, exist_ok=True)
    dataf = get_cluster_members(df)
    dataf.to_csv(sm.params.out, sep="\t")
    shutil.move(sm.params.out, sm.output[0])
    os.rmdir(sm.params.tmpdir)


def main(sm):
    toolbox = {'lulu2tab': lulu2tab}
    toolbox[sm.rule](sm)


if __name__ == '__main__':
    main(snakemake)