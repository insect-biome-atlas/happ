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

## OUTPUT
output <- snakemake@output$tsv

## INPUT
fasta <- snakemake@input$fasta

cat("Reading taxonomy\n")
taxdf <- read.table(snakemake@input$taxonomy, sep="\t", header=TRUE, check.names=FALSE)
# Subset taxonomy to order name
taxonomy <- taxdf[taxdf$Order==order_name,]

countsfile <- snakemake@input$counts
cat("Reading cluster counts\n")
counts <- fread(countsfile, check.names=FALSE, nThread=threads, 
                sep="\t", header=TRUE, data.table=FALSE)
# Subset to clusters defined for order
counts <- counts[counts$cluster %in% taxonomy$cluster,]

# Handle cases where there is only one cluster in the order
if (nrow(taxonomy) == 1) {
     cat("Only one cluster in order\n")
     res <- data.table(matrix(NA, nrow = 1, ncol = 5))
     colnames(res) <- c("cluster", "n_samples", "n_reads", "numt", "reason")
     res$cluster <- taxonomy$cluster
     res$n_samples <- apply(counts[,2:ncol(counts)] > 0, 1, sum)
     res$n_reads <- sum(counts[2:ncol(counts)])
     res$numt <- "FALSE"
     res$reason <- ""
     write.table(res, snakemake@output$tsv, sep="\t", quote=FALSE, row.names=FALSE)
     sink()
     quit(save="no", status=0, runLast=FALSE)
}
res <- data.table(matrix(NA, nrow = 1, ncol = 5))
colnames(res) <- c("cluster", "n_samples", "n_reads", "numt", "reason")

# Add spikeins if not null
if (length(spikeins) > 0) {
     spikes <- taxdf[taxdf$cluster %in% spikeins,]
     taxonomy <- rbind(taxonomy, spikes)
}

res <- combined_filter_neighbors(fasta=fasta, taxonomy=taxonomy, counts=counts, output=output, 
                        spikeins=spikeins, n_closest=n_closest)
sink()