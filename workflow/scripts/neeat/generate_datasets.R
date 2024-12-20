sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)
# Generate dataset partitions
library(data.table)

## RESOURCES
threads <- snakemake@threads

## WILDCARDS
tax <- snakemake@wildcards$tax
rank <- snakemake@wildcards$noise_rank

## INPUT
taxonomy_file <- snakemake@input$taxonomy
counts_file <- snakemake@input$counts

taxonomy <- fread(taxonomy_file, header=TRUE, sep="\t", nThread=threads)
taxonomy <- as.data.frame(taxonomy)
taxonomy <- taxonomy[taxonomy[, "representative"]==1,]
counts <- fread(counts_file, header=TRUE, sep="\t", nThread=threads)

taxdf <- taxonomy[taxonomy[, rank]==tax,]
tax_counts <- counts[counts$cluster %in% taxdf$cluster,]
tax_counts$cluster <- taxdf$ASV[match(tax_counts$cluster,taxdf$cluster)]
write.table(tax_counts, snakemake@output$counts, row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

# Update taxdf to only those clusters that have counts (pos controls are in taxonomy, generating extra names)
# Note that cluster name is now ASV id in counts
taxdf <- taxdf[taxdf$ASV %in% tax_counts$cluster,]
write.table(taxdf, snakemake@output$taxonomy, sep="\t",row.names=FALSE)

sink()