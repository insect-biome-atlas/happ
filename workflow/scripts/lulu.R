sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)

suppressPackageStartupMessages(library(lulu))

## Params
minimum_ratio_type = snakemake@params$minimum_ratio_type
minimum_ratio = snakemake@params$minimum_ratio
minimum_match = snakemake@params$minimum_match
minimum_relative_cooccurence = snakemake@params$minimum_relative_cooccurence

# Read cluster counts
otutab <- read.csv(snakemake@input$otutab, sep='\t',header=TRUE,as.is=TRUE, row.names = 1)
# Check lines in matchlist
matchlist_lines <- length(readLines(snakemake@input$matchlist, n=1))
# if only one cluster in the otutab, or no matches in the matchlist, skip LULU
if (nrow(otutab)==1 | matchlist_lines == 0) {
     cat("Insufficient data (either too few clusters or no pairwise matches). Skipping LULU filtering\n")
     curated_table <- otutab
     otumap <- c()
} else {
     matchlist <- read.table(snakemake@input$matchlist, header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
     curated_result <- lulu(otutab, matchlist, minimum_ratio_type = minimum_ratio_type, minimum_ratio = minimum_ratio, minimum_match = minimum_match, minimum_relative_cooccurence = minimum_relative_cooccurence)
     curated_table <- curated_result$curated_table
     cat("Clusters total: ", nrow(otutab), "\n")
     cat("Clusters retained: ", curated_result$curated_count, "\n")
     cat("Clusters discarded: ", curated_result$discarded_count, "\n")
     cat("Runtime: ", curated_result$runtime, "\n")
     otu_map <- curated_result$otu_map
     # find logfile
     logfile <- list.files(getwd(), pattern="lulu.log")[1]
     file.rename(logfile, snakemake@output$log)
}

# Write details
write.table(otu_map, file=snakemake@output$otu_map, sep="\t", quote=FALSE, row.names=TRUE)

# Write curated table
curated_table <- cbind(rownames(curated_table), curated_table)
colnames(curated_table)[1] <- "ASV"
write.table(curated_table, file=snakemake@output$curated_table, sep="\t", quote=FALSE, row.names=FALSE)

sink()