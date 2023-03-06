#!/bin/bash

module load bioinfo-tools SeqKit biopython cutadapt

python extract_region.py

# Remove gap characters, then remove leading and trailing 'N'
seqkit seq -g finbol.asv-region.fasta | seqkit replace -s -r "" -p "N+$" | seqkit replace -s -r "" -p "^N+" > tmpfile
# Now remove ids still containing non standard DNA chars
seqkit grep -s -r -p "[^ACGTacgt]+" tmpfile | seqkit seq -i | grep ">" | sed 's/>//g' > ids
seqkit grep -v -f ids tmpfile > tmpfile2

cutadapt -m 403 -M 418 -o finbol.asv-region.cleaned.fasta tmpfile2

rm tmpfile tmpfile2 ids
