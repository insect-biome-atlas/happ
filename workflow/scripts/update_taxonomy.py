#/usr/bin/env python

import pandas as pd
import os
from argparse import ArgumentParser
import sys

def assign_lower(items, ranks):
    """
    This function updates the lower ranks so that they are consistent with the updated taxonomy.

    For instance, if the base taxonomy is:
    Kingdom  Phylum     Class   Order                 Family                Genus                   Species
    Animalia Arthropoda Insecta unclassified.Insecta  unclassified.Insecta  unclassified.Insecta    unclassified.Insecta
    and after it has been updated it is:
    Kingdom  Phylum     Class   Order      Family                Genus                   Species
    Animalia Arthropoda Insecta Hemiptera  unclassified.Insecta  unclassified.Insecta    unclassified.Insecta
    Then this function will fix it to:
    Kingdom  Phylum     Class   Order      Family                  Genus                     Species
    Animalia Arthropoda Insecta Hemiptera  unclassified.Hemiptera  unclassified.Hemiptera    unclassified.Hemiptera
    """
    for i, rank in enumerate(ranks):
        if i == 0:
            last_known = items[rank]
            continue
        items[rank] = f"unclassified.{last_known}"
    return items

def main(args):
    # params
    sys.stderr.write(f"Reading updated taxonomy from {args.update}\n")
    for i, update_file in enumerate(args.update):
        if i == 0:
            update = pd.read_csv(update_file, sep="\t", index_col=0, dtype=str)
        else:
            update = pd.concat([update, pd.read_csv(update_file, sep="\t", index_col=0, dtype=str)])
    rename = {}
    if args.rename_file and os.path.exists(args.rename_file):
        sys.stderr.write(f"Mapping queries using {args.rename_file}\n")
        rename_df = pd.read_csv(args.rename_file, sep="\t", index_col=0)
        rename = rename_df.loc[:, args.rename_col].to_dict()
        update = update.assign(ASV=update.index)
        update.rename(index=rename, inplace=True)
        update.index.name=args.rename_col
    sys.stderr.write(f"Reading base taxonomy from {args.base}\n")
    base = pd.read_csv(args.base, sep="\t", index_col=0)
    base.fillna("unclassified", inplace=True)
    # check if update columns should be renamed
    lc_cols = [x for x in base.columns if x.lower() in update.columns]
    if len(lc_cols)>0:
        rename_dict = dict(zip([x.lower() for x in lc_cols], lc_cols))
        update = update.rename(columns=rename_dict)
    # find unclassified entries at update level
    for rank in args.update_ranks:
        unc = base.loc[base[rank].str.startswith("unclassified")].shape[0]
        sys.stderr.write(f"Unclassified at {rank} in base before update: {unc}\n")
    # find clusters that agree on agree_rank level
    update = update.assign(base=base.loc[update.index, args.agree_rank])
    agrees = update.loc[(update[args.agree_rank] == update["base"])&(update[(args.agree_rank)].isin(args.agree_taxa))]
    # Only update entries that are not already classified and are classified in the updated taxonomy
    updated = []
    for rank in args.update_ranks:
        # Get ids that are unclassified in the base taxonomy
        unc_in_base = base.loc[base[rank].str.startswith("unclassified")].index.tolist()
        # Get ids that are classified in the updated taxonomy
        classified_in_update = agrees.loc[~agrees[rank].str.startswith("unclassified")].index.tolist()
        # Ids to update is the intersection between these two
        to_update = list(set(unc_in_base).intersection(classified_in_update))
        sys.stderr.write(f"Entries to update at {rank}: {len(to_update)}\n")
        base.loc[to_update, rank] = agrees.loc[to_update, rank]
        updated+=to_update
    updated = list(set(updated))
    if len(args.downstream_ranks) > 0:
        # Update lower ranks for updated entries
        ranks = [args.update_ranks[0]] + args.downstream_ranks
        sys.stderr.write(f"Updating downstream ranks {' '.join(args.downstream_ranks)}\n")
        base.loc[updated] = base.loc[updated].apply(assign_lower, ranks=ranks, axis=1)
    for rank in args.update_ranks:
        unc = base.loc[base[rank].str.startswith("unclassified")].shape[0]
        sys.stderr.write(f"Unclassified at {rank} in base after update: {unc}\n")
    if not args.output:
        base.to_csv(sys.stdout, sep="\t")
    else:
        base.to_csv(args.output, sep="\t")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-b", "--base", help="Base taxonomy file (to be updated)", required=True)
    parser.add_argument("-u", "--update", help="Updated taxonomy file(s)", nargs="+", required=True)
    parser.add_argument("-a", "--agree_rank", help="Only update entries that agree on <agree_rank>", default="Class")
    parser.add_argument("-t", "--agree_taxa", help="Taxa that must agree on <agree_rank>", nargs="+", default=["Insecta"])
    parser.add_argument("-U", "--update_ranks", help="Ranks to update", nargs="+", default=["Order"])
    parser.add_argument("-d", "--downstream_ranks", help="Ranks downstream of the update ranks for which to update 'unclassified' string", nargs="+")
    parser.add_argument("-r", "--rename_file", help="File used to map ids from update to base")
    parser.add_argument("-c", "--rename_col", help="Column in rename_file to use as id for mapping")
    parser.add_argument("-o", "--output", help="Output file")
    args = parser.parse_args()
    main(args)