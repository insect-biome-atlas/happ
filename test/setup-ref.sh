#!/bin/bash

echo -e "\n\n" >> test/configfile.yml

# get sintax reference
mkdir -p resources/sintax
wget -O resources/sintax/sintax.ref.fasta.gz https://figshare.scilifelab.se/ndownloader/files/38787078
gunzip resources/sintax/sintax.ref.fasta.gz
echo -e "sintax:\n  ref: resources/sintax/sintax.ref.fasta\n" >> test/configfile.yml

# get chesters tree
mkdir -p resources/chesters
wget -O resources/chesters/chesters_new_outgroups_aligned.trim0.9.fasta https://raw.githubusercontent.com/insect-biome-atlas/paper-bioinformatic-methods/refs/heads/main/chesters/data/chesters_tree/chesters_new_outgroups_aligned.trim0.9.fasta
wget -O resources/chesters/chesters_new_outgroups.nwk https://raw.githubusercontent.com/insect-biome-atlas/paper-bioinformatic-methods/refs/heads/main/chesters/data/chesters_tree/chesters_new_outgroups.nwk
wget -O resources/chesters/taxonomy.tsv https://raw.githubusercontent.com/insect-biome-atlas/paper-bioinformatic-methods/refs/heads/main/chesters/data/chesters_tree/taxonomy.tsv

echo -e "epa-ng:\n  tree: resources/chesters/chesters_new_outgroups.nwk\n  ref_taxonomy: resources/chesters/taxonomy.tsv\n  msa: resources/chesters/chesters_new_outgroups_aligned.trim0.9.fasta\n  chunk_size: 100\n" >> test/configfile.yml