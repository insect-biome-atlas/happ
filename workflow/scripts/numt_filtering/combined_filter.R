# Combined algorithm for finding numts and other spurious clusters

sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)

# Libraries needed
library(seqinr)
library(data.table)
# Extra functions
functions <- snakemake@params$functions
c_model <- snakemake@params$codon_model
source(functions)
source(c_model)

## THREADS
threads <- snakemake@threads

## WILDCARDS
order_name <- snakemake@wildcards$order

## PARAMS
spikeins <- readLines(snakemake@input$spikeins)
n_closest <- snakemake@params$n_closest

## INPUT
fasta <- snakemake@input$fasta

cat("Reading taxonomy\n")
taxdf <- read.table(snakemake@input$taxonomy, sep="\t", header=TRUE, check.names=FALSE)
# Subset taxonomy to order name
taxonomy <- taxdf[taxdf$Order==order_name,]

# Add spikeins if not null
if (length(spikeins) > 0) {
     spikes <- taxdf[taxdf$cluster %in% spikeins,]
     taxonomy <- rbind(taxonomy, spikes)
}
countsfile <- snakemake@input$counts
cat("Reading cluster counts\n")
counts <- fread(countsfile, check.names=FALSE, nThread=threads, 
                sep="\t", header=TRUE, data.table=FALSE)
# Subset to clusters defined for order
counts <- counts[counts$cluster %in% taxonomy$cluster,]

## OUTPUT
output <- snakemake@output$tsv
res <- combined_filter_neighbors(fasta=fasta, taxonomy=taxonomy, counts=counts, output=output, 
                        spikeins=spikeins, n_closest=n_closest)
sink()