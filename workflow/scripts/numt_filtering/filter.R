sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)

library(seqinr)

fastafile <- snakemake@input$fasta
taxfile <- snakemake@input$taxonomy
numt_files <- snakemake@input$numt_files
filter_unclassified_rank <- snakemake@params$filter_unclassified_rank

numts <- c()
non_numts <- c()
cat("Reading numt files\n")
for (file in numt_files) {
     df <- read.table(file, sep="\t", header=TRUE, check.names=FALSE)
     numts <- c(numts, df[df$numt == TRUE,"cluster"])
}

cat("Reading taxonomy\n")
taxonomy <- read.table(taxfile, sep="\t", header=TRUE, check.names=FALSE)

cat("Filtering numts\n")
non_numts_df <- taxonomy[!taxonomy$cluster %in% numts,]
if (filter_unclassified_rank%in%colnames(non_numts_df)) {
     cat("Filtering ASVs unclassified at rank:", filter_unclassified_rank, "\n")
     non_numts_df <- non_numts_df[!startsWith(non_numts_df[[filter_unclassified_rank]], "unclassified"),]
}
cat(paste0("Number of non-numt ASVs: ", nrow(non_numts_df), "/", nrow(taxonomy), "\n"))
write.table(non_numts_df, file=snakemake@output$tsv, sep="\t", row.names=FALSE, quote=FALSE)
non_numts_rep_df <- non_numts_df[non_numts_df$representative==1,c("ASV", "cluster")]
cat(paste0("Numer of non-numt clusters: ", nrow(non_numts_rep_df), "/", length(unique(taxonomy$cluster)), "\n"))
non_numts_rep_asvs <- non_numts_rep_df$ASV

seqs <- read.fasta(fastafile, whole.header = FALSE)
non_numts_seqs <- seqs[non_numts_rep_asvs]
seq_names <- paste(non_numts_rep_df$ASV, non_numts_rep_df$cluster, sep=" ")
write.fasta(sequences = non_numts_seqs, names = seq_names, file.out=snakemake@output$fasta)
sink()