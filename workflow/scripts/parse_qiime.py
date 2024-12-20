#!/usr/bin/env python

import pandas as pd
from argparse import ArgumentParser


def extract_taxa(row, taxcol="Taxon", ranks=["kingdom", "phylum", "class", "order", "family", "genus", "species"]):
    rank_translate = {"d__": "domain", "k__": "kingdom", "p__": "phylum", "c__": "class", "o__": "order", "f__": "family", "g__": "genus", "s__": "species", "b__": "bold_bin"}
    data = {}
    items = row[taxcol].split("; ")
    for i, rank in enumerate(ranks):
        try:
            item = items[i][3:]
        except IndexError:
            item = ""
        data[rank] = item
    return pd.Series(data)


def add_unassigned(row):
    _row = row[0:-1]
    if sum(_row!="") == 0:
        return row.replace("", "unclassified")
    last_known = _row[_row!=''].index[-1]
    if last_known == _row.index[-1]:
        return row
    if _row[last_known].startswith("unclassified"):
        uncl_str = _row[last_known]
    else:
        uncl_str = f"unclassified.{_row[_row.index[_row!=''][-1]]}"
    row[row==""] =  uncl_str
    return row


def parse_qiime2(df, taxcol="Taxon", ranks=["kingdom", "phylum", "class", "order", "family", "genus", "species"]):
    """
    Parses the output from QIIME2 with this format:
                                                                                    Taxon  Confidence/Consensus
    Feature ID                                                                                     
    1d17561d16d803f652b0d6bdd671c1cc  k__Bacteria; p__Firmicutes; c__Clostridia; o__...    0.999651
    ac809fd715cced98911f73f1dfb1ffb9  k__Bacteria; p__Proteobacteria; c__Alphaproteo...    0.999861
    ca96b8fd01aa6c508305bacac37da6c2  k__Bacteria; p__Fusobacteria; c__Fusobacteriia...    0.893995
    e1b32002ee69c4e3210a2a257511f96a  k__Bacteria; p__Bacteroidetes; c__Bacteroidia;...    0.992977
    d32d15407d86a2044eead0d72cd3f9e7  k__Bacteria; p__Firmicutes; c__Clostridia; o__...    0.999981

    To a dataframe with separate columns for each rank:
                                    kingdom          phylum                class            order            family          genus species
    Feature ID                                                                                                                               
    1d17561d16d803f652b0d6bdd671c1cc  Bacteria      Firmicutes           Clostridia    Clostridiales   Ruminococcaceae   Oscillospira        
    ac809fd715cced98911f73f1dfb1ffb9  Bacteria  Proteobacteria  Alphaproteobacteria  Rhodobacterales  Rhodobacteraceae     Paracoccus        
    ca96b8fd01aa6c508305bacac37da6c2  Bacteria    Fusobacteria        Fusobacteriia  Fusobacteriales  Fusobacteriaceae  Fusobacterium        
    e1b32002ee69c4e3210a2a257511f96a  Bacteria   Bacteroidetes          Bacteroidia    Bacteroidales                                         
    d32d15407d86a2044eead0d72cd3f9e7  Bacteria      Firmicutes           Clostridia    Clostridiales   Veillonellaceae                       
    """
    parsed = df.apply(extract_taxa, args=(taxcol,ranks,), axis=1).loc[:, ranks].fillna("")
    # add confidence/consensus column
    if "Consensus" in df.columns:
        confcol = "Consensus"
    elif "Confidence" in df.columns:
        confcol = "Confidence"
    parsed = pd.merge(parsed, df.loc[:, confcol], left_index=True, right_index=True)
    return parsed

def main(args):
    df = pd.read_csv(args.input, sep="\t", index_col=0)
    parsed = parse_qiime2(df, taxcol=args.taxcol, ranks=args.ranks)
    parsed = parsed.apply(add_unassigned, axis=1)
    parsed.to_csv(args.output, sep="\t")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("input", help="Input file")
    parser.add_argument("output", help="Output file")
    parser.add_argument("-t", "--taxcol", default="Taxon", help="Taxonomy column (default: Taxon)")
    parser.add_argument("-r", "--ranks", nargs="+", default=["kingdom", "phylum", "class", "order", "family", "genus", "species"], help="Taxonomic ranks (default: kingdom phylum class order family genus species)")
    args = parser.parse_args()
    main(args)