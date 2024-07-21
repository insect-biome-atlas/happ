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
  clustmap$noise_combined <- FALSE
  clustmap$noise_abundance <- FALSE
  rownames(clustmap) <- clustmap$ASV
  clustmap[rownames(filter1[filter1$curated=="merged",]),"noise_combined"] <- TRUE
  clustmap[rownames(filter1[filter1$curated=="merged",]),"noise_abundance"] <- TRUE
  noise_data <- clustmap[, c("cluster","noise_combined", "noise_abundance")]
  rownames(noise_data) <- seq(1, nrow(noise_data))
} else {
  filter1_file <- snakemake@input$filter1
  filter2_file <- snakemake@input$filter2
  filter1 <- read.table(filter1_file, sep="\t", header=TRUE, check.names=FALSE)
  filter2 <- read.table(filter2_file, sep="\t", header=TRUE, check.names=FALSE)
  cat("Combining noise filtering\n")
  # Filter and combine results
  noise_data <- merge(filter1, filter2, by="cluster")
  noise_data <- noise_data[, c("cluster","noise.x", "noise.y")]
  colnames(noise_data) <- c("cluster", "noise_combined", "noise_abundance")
}

noise_data$noise <- ifelse((noise_data$noise_combined == TRUE) | (noise_data$noise_abundance == TRUE), TRUE, FALSE)

## OUTPUT
noise_res_file <- snakemake@output$noise_res
noise_eval_file <- snakemake@output$noise_eval

# Output noise analysis
write.table(x=noise_data, file=noise_res_file, sep="\t", row.names=FALSE, quote=FALSE)

# Execute and output evaluation 
evaluate_res_details(order_name, order_clusters, taxonomy, noise_data, trusted_file, noise_eval_file) 

sink()