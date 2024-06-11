sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)

## PARAMS
functions <- snakemake@params$functions
print(snakemake@params)
source(functions)

## WILDCARDS
order_name <- snakemake@wildcards$order

## PARAMS
trusted_file <- snakemake@params$trusted
lulu <- snakemake@params$lulu

cat("Reading input files\n")
clust_file <- snakemake@input$clust_file
taxfile <- snakemake@input$taxonomy
taxonomy <- read.table(taxfile, sep="\t", check.names=FALSE, header=TRUE)
order_clusters <- read.table(clust_file, sep="\t", header=TRUE, check.names=FALSE)

## INPUT
if (lulu) {
  otu_map <- snakemake@input$otu_map
  # Check if otu_map file is empty, if so then create an empty data frame
  # because we want to keep the single cluster from the order
  if (length(readLines(otu_map, n=2)) == 1) {
    filter1 <- data.frame()
  }else {
    filter1 <- read.table(otu_map, sep="\t", header=TRUE, check.names=FALSE)
  }
  clustmap <- order_clusters[, c("cluster","ASV")]
  clustmap$numt_combined <- FALSE
  clustmap$numt_abundance <- FALSE
  rownames(clustmap) <- clustmap$ASV
  clustmap[rownames(filter1[filter1$curated=="merged",]),"numt_combined"] <- TRUE
  clustmap[rownames(filter1[filter1$curated=="merged",]),"numt_abundance"] <- TRUE
  numt_data <- clustmap[, c("cluster","numt_combined", "numt_abundance")]
  rownames(numt_data) <- seq(1, nrow(numt_data))
} else {
  filter1_file <- snakemake@input$filter1
  filter2_file <- snakemake@input$filter2
  filter1 <- read.table(filter1_file, sep="\t", header=TRUE, check.names=FALSE)
  filter2 <- read.table(filter2_file, sep="\t", header=TRUE, check.names=FALSE)
  cat("Combining numt filtering\n")
  # Filter and combine results
  numt_data <- merge(filter1, filter2, by="cluster")
  numt_data <- numt_data[, c("cluster","numt.x", "numt.y")]
  colnames(numt_data) <- c("cluster", "numt_combined", "numt_abundance")
}

numt_data$numt <- ifelse((numt_data$numt_combined == TRUE) | (numt_data$numt_abundance == TRUE), TRUE, FALSE)

## OUTPUT
numt_res_file <- snakemake@output$numt_res
numt_eval_file <- snakemake@output$numt_eval

# Output numt analysis
write.table(x=numt_data, file=numt_res_file, sep="\t", row.names=FALSE, quote=FALSE)

# Execute and output evaluation 
evaluate_res_details(order_name, order_clusters, taxonomy, numt_data, trusted_file, numt_eval_file) 

sink()