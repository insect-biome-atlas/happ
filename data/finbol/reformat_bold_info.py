#!/usr/bin/env python

import pandas as pd
import sys
from Bio.SeqIO import parse
import re

bold_info = pd.read_csv("bold_info.tsv", sep="\t", index_col=0)
sys.stderr.write(f"# Info for {bold_info.shape[0]} records loaded\n")

bold_info.rename(columns = lambda x: x[0].upper()+x[1:], inplace=True)
bold_info = bold_info.assign(Kingdom=pd.Series(["Animalia"]*bold_info.shape[0], index=bold_info.index))
bold_info = bold_info.loc[:, ["Kingdom","Phylum","Class","Order","Family","Genus","Species","BOLD_bin"]]

# Remove seqs unassigned at Family, Genus, Species or BOLD_bin
asv_taxa = bold_info.loc[~(bold_info.Family!=bold_info.Family)]
asv_taxa = asv_taxa.loc[~(asv_taxa.Genus!=asv_taxa.Genus)]
asv_taxa = asv_taxa.loc[~(asv_taxa.Species!=asv_taxa.Species)]
asv_taxa = asv_taxa.loc[~(asv_taxa.BOLD_bin!=asv_taxa.BOLD_bin)]
sys.stderr.write(f"# {asv_taxa.shape[0]} records remaining after removal of missing values\n")

# Remove Species that do not have "<Genus> <species>" format
sys.stderr.write(f"# Remove records with odd Species name format\n")
def filter_by_regex(genus, species, regex=re.compile("[A-Z][a-z]+ [a-z]+")):
    if genus != species.split(" ")[0]:
        return False
    m = regex.match(species)
    if m==None:
        return False
    if m.group() == species:
        return True
    return False

keep = []
for row in asv_taxa.itertuples():
    if filter_by_regex(row.Genus, row.Species):
        keep.append(row.Index)
sys.stderr.write(f"# Removing {asv_taxa.shape[0]-len(keep)} records with odd Species naming\n")
asv_taxa = asv_taxa.loc[keep]
sys.stderr.write(f"# {asv_taxa.shape[0]} records remaining\n")

# Read and filter sequences
seqs = {} # for storing sequences as keys, record ids as values
recs = {} # for storing records with ids as keys, sequences as values
dupfams = {} # stores number of duplicates per family
sys.stderr.write("# Reading sequences from finbol.asv-region.cleaned.fasta\n")
for record in parse("finbol.asv-region.cleaned.fasta", "fasta"):
    asv = (record.id).split("|")[0]
    if asv in asv_taxa.index:
        recs[asv] = record.seq
        try:
            seqs[str(record.seq)].append(asv)
            fam = asv_taxa.loc[asv, "Family"]
            try:
                dupfams[fam] +=1
            except KeyError:
                dupfams[fam] = 2
        except KeyError:
            seqs[str(record.seq)] = [asv]

sys.stderr.write(f"# {len(recs.keys())} read\n")
asv_taxa = asv_taxa.loc[recs.keys()]
sys.stderr.write(f"# {asv_taxa.shape[0]} records remaining after filtering to records with seqs\n")

dupseqs = len([x for x in seqs.keys() if len(seqs[x]) > 1])
sys.stderr.write(f"# {dupseqs}/{len(seqs.keys())} sequences are non-unique\n")
if dupseqs > 0:
    sys.stderr.write("# Checking ambiguity at Species level\n")
dupasvs = 0 # counts number of 'ASVs' that share these non-unique sequences
ambiguous = [] # list of record ids that are ambiguous at species level
ambig_seqs = 0 # number of sequences that cause ambiguity
with open("ambiguous.tsv", 'w') as fhout:
    fhout.write("seq\tn_records\trecord_ids\tGenus\tSpecies\tBOLD_bin\n")
    for seq, l in seqs.items():
        if len(l) > 1:
            dupasvs+=len(l)
            species = asv_taxa.loc[l, "Species"]
            if len(set(species))>1:
                ambiguous+=l
                ambig_seqs +=1
                n_records = len(l)
                genus = ",".join(asv_taxa.loc[l, "Genus"].unique())
                species = ",".join(asv_taxa.loc[l, "Species"].unique())
                boldbin = ",".join(asv_taxa.loc[l, "BOLD_bin"].unique())
                record_ids = ",".join(l)
                fhout.write(f"{seq}\t{n_records}\t{record_ids}\t{genus}\t{species}\t{boldbin}\n")

sys.stderr.write(f"# non-unique sequences found across {dupasvs} records\n")
if len(ambiguous) > 0:
    sys.stderr.write(f"# {ambig_seqs} sequences have ambiguity at Species level\n")
    sys.stderr.write(f"# {len(ambiguous)} records have ambiguity at Species level\n")
    sys.stderr.write(f"# Removing {len(ambiguous)} records\n")
    asv_taxa = asv_taxa.drop(ambiguous)
    sys.stderr.write(f"# {asv_taxa.shape[0]} records remaining\n")

famsize = asv_taxa.groupby("Family").size()
asv_taxa = asv_taxa.loc[asv_taxa.Family.isin(famsize.loc[famsize>1].index)]
sys.stderr.write(f"# {len(famsize.loc[famsize==1])} Families have only 1 record\n")
sys.stderr.write(f"# {asv_taxa.shape[0]} records remaining after removal of families with just 1 record\n")

sys.stderr.write("# Writing taxa to asv_taxa.tsv\n")
asv_taxa.to_csv("asv_taxa.tsv", sep="\t")
sys.stderr.write("# Writing seqs to asv_seqs.fasta\n")
with open("asv_seqs.fasta", 'w') as fhout:
    for asv, seq in recs.items():
        if asv in asv_taxa.index:
            fhout.write(f">{asv}\n{seq}\n")
