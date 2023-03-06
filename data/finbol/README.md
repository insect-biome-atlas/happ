# FINBOL data

## Obtaining and filtering data

This data was downloaded from BOLDsystems using the link: http://dx.doi.org/10.5883/DS-FINPRO, resulting in `bold_data.txt` and `fasta.fas`.

The `bold_data.txt` file was processed into `bold_info.tsv` like this:

```python
import pandas as pd
data = pd.read_csv("bold_data.txt", sep="\t", index_col=0, header=0, dtype=str)
cols = ["bin_uri", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name"]
df = data.loc[:, cols]
df.columns = ["BOLD_bin", "phylum", "class", "order", "family", "genus", "species"]
df.to_csv("bold_info.tsv", sep="\t")
```

The `fasta.fas` file was filtered using the `filter_COI.py` script resulting in `finbol.COI.fasta`.

## Alignments

Sequences in the `finbol.COI.fasta` file were aligned using `clustalo` using the `msa.sh` batch script:

```bash
module load bioinfo-tools clustalo/1.2.4

fasta="finbol.COI.fasta"
out="finbol.msa.nohmm.fasta"

clustalo -i $fasta -o $out --dealign -t DNA --threads 20 -l clustalo.nohmm.log -v
```

Using the MSAs from clustal omega I checked sequence LEFIA236-10 in FINBOL assigned as Selenia dentaria,
then looked up an ASV with the same species assigned in the MAD data (6da411eda0d9c8e9ba72acae3981b062)
and searched for the starting/ending sequence of that ASV in the alignment (opened with Aliview).
That identified position 284 - 750 in the alignment.

Using FIDIP3710-13 (Brachypeza bisignata) and 14520346d04ff4e16c4aa4f06a156bd1 in the MBs_merged data
gave the same positions.

## Extracting the ASV region
The ASV region corresponding to positions 284-750 in the alignment was extracted with `extract_region.py` 
yielding `finbol.asv-region.fasta`. This was followed by removal of gap-characters and leading and trailing
'N' and non standard DNA characters, as in the coidb-tool.

```bash
module load bioinfo-tools SeqKit biopython cutadapt

python extract_region.py

# Remove gap characters, then remove leading and trailing 'N'
seqkit seq -g finbol.asv-region.fasta | seqkit replace -s -r "" -p "N+$" | seqkit replace -s -r "" -p "^N+" > tmpfile
# Now remove ids still containing non standard DNA chars
seqkit grep -s -r -p "[^ACGTacgt]+" tmpfile | seqkit seq -i | grep ">" | sed 's/>//g' > ids
seqkit grep -v -f ids tmpfile > tmpfile2

cutadapt -m 403 -M 418 -o finbol.asv-region.cleaned.fasta tmpfile2

rm tmpfile tmpfile2 ids
```

This resulted in `finbol.asv-region.cleaned.fasta`

## Taxonomy
Next came the issues with taxonomy. I first found that 741 (of 38267) entries in the the Finbol database did not have
Species or BOLD_bin assignments (680 lacked species, 69 lacked BOLD_bin and 8 lacked both). These were all excluded.
Further 35 records were removed because they were the only record within their Family.

I also found that some sequences in the fasta file were identical but came from records that differed in taxonomy, even
between species. After some discussion we concluded that the reason for this is that the ASV fragment we used is not
long enough to differentiate between sequences. Therefore we decided to exclude records corresponding to ASV sequences 
with multiple species. 

The `bold_info.tsv` and `finbol.asv-region.cleaned.fasta` files were processed using the `reformat_bold_info.py` script
located at `/crex/proj/snic2020-16-248/nobackup/ASV-clustering/data/finbol`. The code for that script is also at the 
bottom of this README.

## Appendix

### Scripts

`extract_region.py`

```python
#!/usr/bin/env python
###################################
from Bio.SeqIO import parse

infile = "finbol.msa.nohmm.fasta"
outfile = "finbol.asv-region.fasta"
with open(outfile, 'w') as fhout:
    for record in parse(infile, 'fasta'):
        seq = record.seq[283:750]
        if seq.startswith("-"):
            continue
        fhout.write(f">{record.description}\n{seq}\n")
```

`filter_COI.py`

```python
#!/usr/bin/env python
###################################
from Bio.SeqIO import parse

with open("finbol.COI.fasta", 'w') as fhout:
    for record in parse("fasta.fas", 'fasta'):
        desc = record.description
        gene = desc.split("|")[2]
        if gene == "COI-5P":
            fhout.write(f">{desc}\n{record.seq}\n")
```

`reformat_bold_info.py`

```python
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
```

**Log output from reformat_bold_info.py**:

```
# Info for 38275 records loaded
# 37534 records remaining after removal of missing values
# Remove records with odd Species name format
# Removing 1954 records with odd Species naming
# 35580 records remaining
# Reading sequences from finbol.asv-region.cleaned.fasta
# 32401 read
# 32401 records remaining after filtering to records with seqs
# 6434/18680 sequences are non-unique
# Checking ambiguity at Species level
# non-unique sequences found across 20155 records
# 349 sequences have ambiguity at Species level
# 1439 records have ambiguity at Species level
# Removing 1439 records
# 30962 records remaining
# 33 Families have only 1 record
# 30929 records remaining after removal of families with just 1 record
# Writing taxa to asv_taxa.tsv
# Writing seqs to asv_seqs.fasta
```


