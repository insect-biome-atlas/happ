sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)

# Libraries needed
library(data.table)

## Extra functions
functions <- snakemake@params$functions
source(functions)

## Params
abundance_threshold <- snakemake@params$abundance_threshold

## WILDCARDS
order_name <- snakemake@wildcards$order

threads <- snakemake@threads

## INPUT
taxonomy <- read.table(snakemake@input$taxonomy, sep="\t", header=TRUE, check.names=FALSE)
# Subset taxonomy to order name
taxonomy <- taxonomy[taxonomy$Order==order_name,]
countsfile <- snakemake@input$counts
cat("Reading cluster counts\n")
counts <- fread(countsfile, check.names=FALSE, nThread=threads, 
                sep="\t", header=TRUE, data.table=FALSE)
# Subset to clusters defined for order
counts <- counts[counts$cluster %in% taxonomy$cluster,]

## OUTPUT
output <- snakemake@output$tsv
res <- abundance_filter(order_name=order_name, taxonomy=taxonomy, counts=counts, output=output, save_to_file=TRUE, threshold=abundance_threshold)
sink()
