# Script for computing evolutionary distance lists (including pdist)
library(seqinr)
source("codon_model.R")
source("get_data_fxns.R")

matchlist <- read.delim("../bigdata/vsearch_matchlist.tsv",header=FALSE)

T <- get_se_cluster_taxonomy()
T <- T[T$representative==1 & !grepl("_X",T$Order) & !grepl("unclassified",T$Order),]

orders <- unique(T$Order)

for (ord in orders ) {

    cat("Processing ",ord,"\n",sep="")
    cat ("0%                                                100%\n")
    cat ("|")

    ord_tax <- T[T$Order==ord,]
    ord_asvs <- ord_tax$ASV
    seqs <- read.alignment(paste0("../alignments/",ord,"_aligned.fasta"),format="fasta")

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
            
            if ((s1_name %in% ord_asvs) && (s2_name %in% ord_asvs)) {
                seq1 <- seqs$seq[[which(seqs$nam==s1_name)]]
                seq2 <- seqs$seq[[which(seqs$nam==s2_name)]]

                evodistlist$dadn[i] <- dadn_dist(seq1, seq2)
                evodistlist$wdadn[i] <- wdadn_dist(seq1, seq2)
            }
        }
    }

    cat("|\n")

    write.table(evodistlist, paste0("../data/",ord,"_evodistlist.tsv"), row.names=FALSE, sep="\t")
}

