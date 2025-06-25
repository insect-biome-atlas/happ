This folder contains a small test dataset which can be used to test the standalone neeat pipeline. A set of 100 clusters/OTUs were selected from the full run with HAPP on the [Swedish IBA data](https://doi.org/10.1038/s41597-025-05151-0). 


## File descriptions

- `cluster_counts.tsv`: A tab-separated file with counts for 100 clusters/OTUs (rows) in 102 samples (columns)
- `cluster_reps.fasta`: A fasta file with sequences of the 100 cluster/OTU representatives
- `cluster_taxonomy.tsv`: A tab-separated file containing taxonomic assignments for all sequences as well as their cluster/OTU assignments. Rows with the `representative` column set to 1 indicate sequences which are representatives for their cluster/OTU.
- `metadata.tsv`: A tab-separated file with sample ids in the first column and sample types in the second column. This file is supplied as an example of how sample metadata can be used in HAPP to only use counts for a certain sample type.