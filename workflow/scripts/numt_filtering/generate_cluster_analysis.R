sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)

## LIBRARIES
library(data.table)
functions <- snakemake@params$functions
source(functions)

## INPUT
taxfile <- snakemake@input$taxonomy
countsfile <- snakemake@input$counts

## PARAMS
consensus_taxfile <- snakemake@params$consensus_taxonomy

## OUTPUT
outfile <- snakemake@output$tsv

## WILDCARDS
order_name <- snakemake@wildcards$order

threads <- snakemake@threads

complementing_data(order_name, taxfile, consensus_taxfile, countsfile, outfile, threads)

sink()