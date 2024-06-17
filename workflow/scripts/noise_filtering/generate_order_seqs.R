#Generate unaligned sequence files for orders

sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)

# Libraries needed
library(data.table)
library(seqinr)

# Read in data on rep asvs so we can translate to clusters
cat("Reading nochimera cluster info\n")
taxonomy <- read.table(snakemake@input$taxonomy,header=TRUE,sep="\t")

# Read in fasta alignment 
cat("Reading alignment\n")
seqs <- read.fasta(snakemake@input$fasta, whole.header = FALSE)

# Get list of orders
#orders_insecta <- unique(taxonomy$Order[taxonomy$Class=="Insecta"])
#orders_collembola <- unique(taxonomy$Order[taxonomy$Class=="Collembola"])
#orders <- c(orders_insecta, orders_collembola)
order_name <- snakemake@wildcards$order
rep_asvs <- taxonomy$ASV[taxonomy$Order==order_name & taxonomy$representative==1]
rep_seqs <- seqs[rep_asvs]
write.fasta(sequences = rep_seqs, names = names(rep_seqs), file.out=snakemake@output[[1]])

sink()