#!/usr/bin/env python

from argparse import ArgumentParser
from Bio import Phylo


def main(args):
    ranks = args.ranks
    tree = Phylo.read(args.tree, format="newick")
    with open(args.outfile, 'w') as fhout:
        for t in tree.get_terminals():
            name = t.name
            taxnames = ";".join(name.split("_")[0:len(ranks)-1]+["_".join(name.split("_")[len(ranks)-2:])])
            fhout.write(f"{name}\t{taxnames}\n")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("tree", type=str,
                        help="Input tree")
    parser.add_argument("outfile", type=str,
                        help="Output taxon-file")
    parser.add_argument("--ranks", nargs="*", default=["order","family","genus","species"],
                        help="Ranks to extract")
    args = parser.parse_args()
    main(args)