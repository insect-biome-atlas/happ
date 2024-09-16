library(ape)

data_path <- "~/dev/figshare-repos/iba/processed_data/SE.v2/"

seqs <- read.FASTA(paste0(data_path,"cluster_reps.fasta"))

# The names will be "<asv> <cluster_name>". Convert this to
# "<asv>".
x <- strsplit(names(seqs), " ", fixed=TRUE)
y <- character(length(x))
for (i in 1:length(x)) {
    y[i] <- unlist(x[i])[1]
}
names(seqs) <- y

seqs_aa <- trans(seqs, code=5, codonstart=2)

write.FASTA(seqs_aa, "../alignments/cluster_reps_aa.fasta")

# Write order-specific files
T <- read.delim("~/dev/figshare-repos/iba/processed_data/SE.v2/cluster_taxonomy.tsv",header=TRUE,sep="\t")
T <- T[T$representative==1 & !grepl("_X", T$Order) & !grepl("unclassified", T$Order),]

orders <- unique(T$Order)
for (ord in orders) {

    ord_asvs <- T$ASV[T$Order==ord]
    ord_seqs <- seqs[which(names(seqs) %in% ord_asvs)]
    ord_seqs_aa <- seqs_aa[which(names(seqs_aa) %in% ord_asvs)]

    write.FASTA(ord_seqs, file=paste0("../alignments/",ord,".fasta"))
    write.FASTA(ord_seqs_aa, file=paste0("../alignments/",ord,"_aa.fasta"))
}

write.table(data.frame(orders), file="../alignments/orders.tsv",sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

