sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)

## PARAMS
functions <- snakemake@params$functions
print(snakemake@params)
source(functions)

## WILDCARDS
order_name <- snakemake@wildcards$order

## INPUT
filter1_file <- snakemake@input$filter1
filter2_file <- snakemake@input$filter2
clust_file <- snakemake@input$clust_file
taxfile <- snakemake@input$taxonomy

## PARAMS
trusted_file <- snakemake@params$trusted

## OUTPUT
numt_res_file <- snakemake@output$numt_res
numt_eval_file <- snakemake@output$numt_eval

cat("Reading input files\n")
filter1 <- read.table(filter1_file, sep="\t", header=TRUE, check.names=FALSE)
filter2 <- read.table(filter2_file, sep="\t", header=TRUE, check.names=FALSE)
order_clusters <- read.table(clust_file, sep="\t", header=TRUE, check.names=FALSE)
taxonomy <- read.table(taxfile, sep="\t", check.names=FALSE, header=TRUE)

cat("Combining numt filtering\n")
# Filter and combine results
numt_data <- merge(filter1, filter2, by="cluster")
numt_data <- numt_data[, c("cluster","numt.x", "numt.y")]
colnames(numt_data) <- c("cluster", "numt_combined", "numt_abundance")
numt_data$numt <- ifelse((numt_data$numt_combined == TRUE) | (numt_data$numt_abundance == TRUE), TRUE, FALSE)

# Output numt analysis
write.table(x=numt_data, file=numt_res_file, sep="\t", row.names=FALSE, quote=FALSE)

# Execute and output evaluation 
evaluate_res_details(order_name, order_clusters, taxonomy, numt_data, trusted_file, numt_eval_file) 

sink()