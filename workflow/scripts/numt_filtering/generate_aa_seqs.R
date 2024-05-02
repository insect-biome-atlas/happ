#Generate unaligned sequence files for each insect order

sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)

# Libraries needed
library(ape)
library(data.table)

# Read in data on rep asvs so we can translate to clusters
#cat("Reading nochimera cluster info\n")
#taxonomy <- read.table(snakemake@input$taxonomy,header=TRUE,sep="\t")

# Read in fasta alignment 
cat("Reading sequences\n")
seqs <- read.FASTA(snakemake@input$fasta)

# Get list of orders
#orders_insecta <- unique(taxonomy$Order[taxonomy$Class=="Insecta"])
#orders_collembola <- unique(taxonomy$Order[taxonomy$Class=="Collembola"])
#orders <- c(orders_insecta, orders_collembola)

order_name <- snakemake@wildcards$order
seqs <- read.FASTA(file=snakemake@input$fasta)
aa_seqs <- trans(seqs, code=5)
write.FASTA(aa_seqs, file=snakemake@output[[1]])

sink()