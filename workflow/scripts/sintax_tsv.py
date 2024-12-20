#!/usr/bin/env python3
##### Program description #######
#
# Title: Parse SINTAX output
#
# Author(s): Lokeshwaran Manoharan
#
#
#
# Description: Parse the SINTAX output and create a tsv file with taxonomy
#  
# List of subroutines: 
#
#
#
# Overall procedure: using functions and dictionaries in python
#
# Usage: sintax_tsv.py -i <sintax output> -o <sintax tsv output>
##################################

import re
import sys
import argparse

usage = '''This program takes the output of SINTAX and creates readable taxonomy file'''

parser = argparse.ArgumentParser(description=usage)


parser.add_argument(
    '-i', '--infile',
    dest='infile',
    metavar='INFILE',
    type=argparse.FileType('r'),
    help='SINTAX output',
    required=True
    )

parser.add_argument(
    '-o', '--outfile',
    dest='outfile',
    metavar='OUTFILE',
    type=argparse.FileType('w'),
    help='Final Taxonomy in tsv format',
    default=sys.stdout
    )

parser.add_argument("-r", "--ranks", nargs="+", help="List of ranks to include in the output", 
                    default=["kingdom", "phylum", "class", "order", "family", "genus", "species"])

args = parser.parse_args()

def complete_list(some_list, target_len):
    return some_list[:target_len] + ['unclassified.' + some_list[-1]]*(target_len - len(some_list))

p0 = re.compile('>')
p1 = re.compile(' ')
p2 = re.compile(',')
p3 = re.compile('\t')
p4 = re.compile(';')

header = ["ASV"]+args.ranks

print("\t".join(header), sep = '', file = args.outfile)

for line in args.infile:
    line = line.rstrip('\n')
    tmp_list = re.split(p3, line)
    if tmp_list[3] != '':
        annot = re.sub('[a-z]:', '', tmp_list[3])
        tmp_list1 = re.split(p2, annot)
        if len(tmp_list1) != len(args.ranks):
            tmp_list1 = complete_list(tmp_list1, len(args.ranks))
            print(tmp_list[0], '\t'.join(tmp_list1), sep = '\t', file = args.outfile)
        else:
            print(tmp_list[0], '\t'.join(tmp_list1), sep = '\t', file = args.outfile)
    else:
        print(tmp_list[0], '\t'.join(['unclassified'] * len(args.ranks)), sep = '\t', file = args.outfile)

args.outfile.close()
args.infile.close()