sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)
library(ape)

## Params
codon_table <- snakemake@params$codon_table
codon_start <- snakemake@params$codon_start

seqs <- read.FASTA(snakemake@input$fasta)

# The names will be "<asv> <cluster_name>". Convert this to
# "<asv>".
x <- strsplit(names(seqs), " ", fixed=TRUE)
y <- character(length(x))
for (i in 1:length(x)) {
    y[i] <- unlist(x[i])[1]
}
names(seqs) <- y

seqs_aa <- trans(seqs, code=codon_table, codonstart=codon_start)

write.FASTA(seqs_aa, snakemake@output$faa)

sink()