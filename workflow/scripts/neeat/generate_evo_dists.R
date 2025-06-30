sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)

# Script for computing evolutionary distance lists (including pdist)

library(seqinr)
source(snakemake@params$codon_model)

matchlist_file <- snakemake@input$matchlist
taxonomy_file <- snakemake@input$taxonomy
fasta_file <- snakemake@input$fasta
tax <- snakemake@wildcards$tax

# Check lines in matchlist
matchlist_lines <- length(readLines(matchlist_file, n=1))
if (matchlist_lines == 0) {
     cat("Insufficient data (no pairwise matches)\n")
     file.create(snakemake@output$tsv)
     quit(save="no", status=0)
}
matchlist <- read.delim(matchlist_file,header=FALSE)

T <- read.delim(taxonomy_file)
T <- T[T$representative==1,]
rownames(T) <- T$cluster

tax_asvs <- T$ASV
seqs <- read.alignment(fasta_file,format="fasta")

evodistlist <- matchlist[(matchlist$V1 %in% seqs$nam) & (matchlist$V2 %in% seqs$nam),]
colnames(evodistlist) <- c("asv1","asv2","idty")
evodistlist$pdist <- 1.0 - (evodistlist$idty / 100.0)
evodistlist$dadn <- rep(NA,nrow(evodistlist))
evodistlist$wdadn <- rep(NA,nrow(evodistlist))
print_interval <- ceiling(nrow(evodistlist) / 50)

if (nrow(evodistlist>0)) {
    for (i in 1:nrow(evodistlist)) {
        if (i%%print_interval==0)
            cat("-")
        s1_name <- evodistlist$asv1[i]
        s2_name <- evodistlist$asv2[i]
        if ((s1_name %in% tax_asvs) && (s2_name %in% tax_asvs)) {
            seq1 <- seqs$seq[[which(seqs$nam==s1_name)]]
            seq2 <- seqs$seq[[which(seqs$nam==s2_name)]]
            evodistlist$dadn[i] <- dadn_dist(seq1, seq2)
            evodistlist$wdadn[i] <- wdadn_dist(seq1, seq2)
        }
    }
    cat("|\n")
}

write.table(evodistlist, snakemake@output$tsv, row.names=FALSE, sep="\t")

sink()