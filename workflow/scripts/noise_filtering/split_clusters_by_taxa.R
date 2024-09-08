sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)

## LIBRARIES
library(data.table)
functions <- snakemake@params$functions
source(functions)

## INPUT
taxfile <- snakemake@input$taxonomy
countsfile <- snakemake@input$counts
consensus_taxfile <- snakemake@input$consensus_taxonomy
fastafile <- snakemake@input$fasta
## OUTPUT
outdir <- snakemake@output[[1]]

## WILDCARDS
order_name <- snakemake@wildcards$order

threads <- snakemake@threads

# Read in cluster information and create core of the data
# on the order clusters
cat("Reading cluster info\n")
cluster_taxonomy <- read.table(taxfile,header=TRUE,sep="\t") # taxonomy
colnames(cluster_taxonomy) <- tolower(colnames(cluster_taxonomy))
orders <- unique(cluster_taxonomy$order)

# Read in the consensus taxonomy
cat("Reading consensus taxonomy\n")
consensus_taxonomy <- read.table(consensus_taxfile,header=TRUE,sep="\t")
colnames(consensus_taxonomy) <- tolower(colnames(consensus_taxonomy))

# Get total number of reads for the clusters
cat("Reading cluster counts\n")
cluster_counts <- fread(countsfile, check.names=FALSE, nThread=threads, 
          sep="\t", header=TRUE, data.table=FALSE, )

# Read in fasta
cat("Reading fasta\n")
seqs <- read.fasta(fastafile, whole.header = FALSE)

for (order_name in unique(cluster_taxonomy$order)) {
     outfile <- paste0(outdir, "/", order_name, "_cluster_analysis.tsv")
     complementing_data(order_name, cluster_taxonomy, consensus_taxonomy, cluster_counts, outfile)
     rep_asvs <- cluster_taxonomy$asv[cluster_taxonomy$order==order_name & cluster_taxonomy$representative==1]
     rep_seqs <- seqs[rep_asvs]
     outfile <- paste0(outdir, "/", order_name, "_seqs.fasta")
}

sink()