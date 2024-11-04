#!/usr/bin/env python

from argparse import ArgumentParser
from Bio import Phylo
import sys

def rename(tree, mapfile):
    m = {}
    with open(mapfile, 'r') as fhin:
        for line in fhin:
            old, new = line.rstrip().split("\t")
            m[old] = new
    for terminal in tree.get_terminals():
        try:
            terminal.name = m[terminal.name]
        except KeyError:
            sys.exit(f"Could not find {terminal.name} in mapfile\n")
    return tree


def main(args):
    tree = Phylo.read(file=args.infile, format="nexus")
    if args.mapfile:
        tree = rename(tree, args.mapfile)
    Phylo.write(trees=tree, file=args.outfile, format="newick")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("infile", type=str,
                        help="Input tree in nexus format")
    parser.add_argument("outfile", type=str,
                        help="Output tree in newick format")
    parser.add_argument("--mapfile", type=str,
                        help="Mapping file used to rename terminals")
    args = parser.parse_args()
    main(args)