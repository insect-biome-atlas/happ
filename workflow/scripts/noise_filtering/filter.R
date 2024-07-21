sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)

library(seqinr)

fastafile <- snakemake@input$fasta
taxfile <- snakemake@input$taxonomy
noise_files <- snakemake@input$noise_files
filter_unclassified_rank <- snakemake@params$filter_unclassified_rank

noise <- c()
non_noise <- c()
cat("Reading noise files\n")
for (file in noise_files) {
     df <- read.table(file, sep="\t", header=TRUE, check.names=FALSE)
     noise <- c(noise, df[df$noise == TRUE,"cluster"])
}

cat("Reading taxonomy\n")
taxonomy <- read.table(taxfile, sep="\t", header=TRUE, check.names=FALSE)

cat("Filtering noise\n")
non_noise_df <- taxonomy[!taxonomy$cluster %in% noise,]
if (filter_unclassified_rank%in%colnames(non_noise_df)) {
     cat("Filtering ASVs unclassified at rank:", filter_unclassified_rank, "\n")
     non_noise_df <- non_noise_df[!startsWith(non_noise_df[[filter_unclassified_rank]], "unclassified"),]
}
cat(paste0("Number of non-noise ASVs: ", nrow(non_noise_df), "/", nrow(taxonomy), "\n"))
write.table(non_noise_df, file=snakemake@output$tsv, sep="\t", row.names=FALSE, quote=FALSE)
non_noise_rep_df <- non_noise_df[non_noise_df$representative==1,c("ASV", "cluster")]
cat(paste0("Numer of non-noise clusters: ", nrow(non_noise_rep_df), "/", length(unique(taxonomy$cluster)), "\n"))
non_noise_rep_asvs <- non_noise_rep_df$ASV

seqs <- read.fasta(fastafile, whole.header = FALSE)
non_noise_seqs <- seqs[non_noise_rep_asvs]
seq_names <- paste(non_noise_rep_df$ASV, non_noise_rep_df$cluster, sep=" ")
write.fasta(sequences = non_noise_seqs, names = seq_names, file.out=snakemake@output$fasta)
sink()