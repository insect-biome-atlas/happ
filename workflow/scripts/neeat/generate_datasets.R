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
cal_counts_file <- snakemake@input$cal_counts
tot_prop_counts_file <- snakemake@input$tot_prop_counts
sample_prop_counts_file <- snakemake@input$sample_prop_counts

taxonomy <- fread(taxonomy_file, header=TRUE, sep="\t", nThread=threads)
taxonomy <- as.data.frame(taxonomy)
taxonomy <- taxonomy[taxonomy[, "representative"]==1,]
counts <- fread(counts_file, header=TRUE, sep="\t", nThread=threads)
cal_counts <- fread(cal_counts_file, header=TRUE, sep="\t", nThread=threads)
tprop_counts <- fread(tot_prop_counts_file, header=TRUE, sep="\t", nThread=threads)
sprop_counts <- fread(sample_prop_counts_file, header=TRUE, sep="\t", nThread=threads)

taxdf <- taxonomy[taxonomy[, rank]==tax,]
tax_counts <- counts[counts$cluster %in% taxdf$cluster,]
tax_counts$cluster <- taxdf$ASV[match(tax_counts$cluster,taxdf$cluster)]
write.table(tax_counts, snakemake@output$counts, row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

# Update taxdf to only those clusters that have counts (pos controls are in taxonomy, generating extra names)
# Note that cluster name is now ASV id in counts
taxdf <- taxdf[taxdf$ASV %in% tax_counts$cluster,]
write.table(taxdf, snakemake@output$taxonomy, sep="\t",row.names=FALSE)

tax_counts <- cal_counts[cal_counts$cluster %in% taxdf$cluster,]
tax_counts$cluster <- taxdf$ASV[match(tax_counts$cluster,taxdf$cluster)]
write.table(tax_counts, snakemake@output$cal_counts, row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

tax_counts <- tprop_counts[tprop_counts$cluster %in% taxdf$cluster,]
tax_counts$cluster <- taxdf$ASV[match(tax_counts$cluster,taxdf$cluster)]
write.table(tax_counts, snakemake@output$tot_prop_counts, row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

tax_counts <- sprop_counts[sprop_counts$cluster %in% taxdf$cluster,]
tax_counts$cluster <- tax_counts$ASV[match(tax_counts$cluster,taxdf$cluster)]
write.table(tax_counts, snakemake@output$sample_prop_counts, row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

sink()