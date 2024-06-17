#Generate unaligned sequence files for each insect order

sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)

# Libraries needed
library(seqinr)
library(data.table)

# Read in data on rep asvs so we can translate to clusters
#cat("Reading nochimera cluster info\n")
#taxonomy <- read.table(snakemake@input$taxonomy,header=TRUE,sep="\t")

# Get list of orders
#orders_insecta <- unique(taxonomy$Order[taxonomy$Class=="Insecta"])
#orders_collembola <- unique(taxonomy$Order[taxonomy$Class=="Collembola"])
#orders <- c(orders_insecta, orders_collembola)

order_name <- snakemake@wildcards$order
cat("Reading sequences\n")
seqs <- read.alignment(file=snakemake@input$fasta, format="fasta")

cat("Trimming sequences\n")
for (i in 1:length(seqs$seq)) {
    seqs$seq[[i]] <- substring(seqs$seq[[i]],first=2)
  }

write.fasta(seqs$seq, seqs$nam, file.out=snakemake@output[[1]], nbchar=500)

sink()