library(data.table)

# Help function for summarizing row reads
row_reads <- function(dt) {

    x <- numeric(nrow(dt))
    for (i in 1:nrow(dt)) {
        if (i %% 10000 == 0)
            cat("Processed",i,"rows (",round(100*i/nrow(dt),digits=1),"%)\n")
        x[i] <- sum(as.numeric(dt[i,-1]))
    }

    x
}

# Matching total number of reads to the actual clusters
match_reads <- function(dt,dt_reads,target_clusters) {

    idx <- match(target_clusters, unlist(dt[,1]))
    dt_reads[idx]
}

# Help function to obtain reads and counts
num_resolved <- function(dt, rank, dt_reads) {

    names <- unlist(dt[,..rank])
    x <- grepl("unclassified",names) | grepl("_X",names) | names==""

    list(count=sum(!x), reads=sum(dt_reads[!x]))
}

# Help function to add results for a method, country and count type
add_results <- function(df, method, country, count_type, A, A_reads) {

    prefix <- list(method=method, country=country, count_type=count_type)

    for (rank in c("Phylum","Class","Order","Family","Genus","Species")) {
        
        res <- num_resolved(A, rank, A_reads)

        if (rank=="Phylum") {
            df <- rbind(df, c(prefix,
                              rank=rank,
                              count=res$count,
                              reads=res$reads,
                              prop_count=1.0,
                              prop_reads=1.0))
            phylum_count <- res$count
            phylum_reads <- res$reads
        } else {
            df <- rbind(df, c(prefix,
                              rank=rank,
                              count=res$count,
                              reads=res$reads,
                              prop_count=res$count/phylum_count,
                              prop_reads=res$reads/phylum_reads))
        }
    }

    df
}

annotation_res <- data.frame(list(method=character(),
                                  country=character(),
                                  count_type=character(),
                                  rank=character(),
                                  count=numeric(),
                                  reads=numeric(),
                                  prop_count=numeric(),
                                  prop_reads=numeric()))

# Analyze Swedish annotations at ASV level
cat("Reading ASV counts for SE\n")
C <- fread("~/dev/figshare-repos/iba/raw_data/CO1_asv_counts_SE.tsv")
cat("Computing reads for SE ASVs\n")
C_reads <- row_reads(C)
cat("done\n")

S1<-fread("../bigdata/sintax_taxonomy_SE.tsv")
S1_reads <- match_reads(C,C_reads,S1$ASV)
annotation_res <- add_results(annotation_res, "sintax", "Sweden", "ASVs", S1, S1_reads)

V1<-fread("../bigdata/vsearch_taxonomy_SE.tsv")
colnames(V1) <- c("ASV","Kingdom","Phylum","Class","Order","Family","Genus","Species","Consensus")
V1_reads <- match_reads(C,C_reads,V1$ASV)
annotation_res <- add_results(annotation_res, "vsearch", "Sweden", "ASVs", V1, V1_reads)

E1<-fread("../bigdata/epang_baseball_taxonomy_SE.tsv")
E1_reads <- match_reads(C,C_reads,E1$ASV)
annotation_res <- add_results(annotation_res, "epang", "Sweden", "ASVs", E1, E1_reads)

EE1<-fread("../bigdata/epang_baseball_taxonomy_chesters_expanded_SE.tsv")
EE1_reads <- match_reads(C,C_reads,EE1$ASV)
annotation_res <- add_results(annotation_res, "epang_exp", "Sweden", "ASVs", EE1, EE1_reads)

rm(C)

# Analyze Malagasy annotations at ASV level
cat("Reading ASV counts for MG\n")
C <- fread("~/dev/figshare-repos/iba/raw_data/CO1_asv_counts_MG.tsv")   # Sample labels are irrelevant
cat("Computing reads for MG ASVs\n")
C_reads <- row_reads(C)
cat("done\n")

S2<-fread("../bigdata/sintax_taxonomy_MG.tsv")
S2_reads <- match_reads(C,C_reads,S2$ASV)
annotation_res <- add_results(annotation_res, "sintax", "Madagascar", "ASVs", S2, S2_reads)

V2<-fread("../bigdata/vsearch_taxonomy_MG.tsv")
colnames(V2) <- c("ASV","Kingdom","Phylum","Class","Order","Family","Genus","Species","Consensus")
V2_reads <- match_reads(C,C_reads,V2$ASV)
annotation_res <- add_results(annotation_res, "vsearch", "Madagascar", "ASVs", V2, V2_reads)

E2<-fread("../bigdata/epang_baseball_taxonomy_MG.tsv")
E2_reads <- match_reads(C,C_reads,E2$ASV)
annotation_res <- add_results(annotation_res, "epang", "Madagascar", "ASVs", E2, E2_reads)

EE2<-fread("../bigdata/epang_baseball_taxonomy_chesters_expanded_MG.tsv")
EE2_reads <- match_reads(C,C_reads,EE2$ASV)
annotation_res <- add_results(annotation_res, "epang_exp", "Madagascar", "ASVs", EE2, EE2_reads)

rm(C)

# Analyze Swedish annotations at cluster level
C <- fread("~/dev/figshare-repos/iba/processed_data/SE.v2/cluster_counts.tsv")
C_reads <- row_reads(C)
T1 <- read.delim("~/dev/figshare-repos/iba/processed_data/SE.v2/cluster_taxonomy.tsv")
T1 <- T1[T1$representative==1,]

S1 <- S1[S1$ASV %in% T1$ASV,]
S1_reads <- match_reads(C,C_reads,T1$cluster[match(S1$ASV,T1$ASV)])
annotation_res <- add_results(annotation_res, "sintax", "Sweden", "clusters", S1, S1_reads)

V1 <- V1[V1$ASV %in% T1$ASV,]
V1_reads <- match_reads(C,C_reads,T1$cluster[match(V1$ASV,T1$ASV)])
annotation_res <- add_results(annotation_res, "vsearch", "Sweden", "clusters", V1, V1_reads)

E1 <- E1[E1$ASV %in% T1$ASV,]
E1_reads <- match_reads(C,C_reads,T1$cluster[match(E1$ASV,T1$ASV)])
annotation_res <- add_results(annotation_res, "epang", "Sweden", "clusters", E1, E1_reads)

EE1 <- EE1[EE1$ASV %in% T1$ASV,]
E1_reads <- match_reads(C,C_reads,T1$cluster[match(EE1$ASV,T1$ASV)])
annotation_res <- add_results(annotation_res, "epang_exp", "Sweden", "clusters", EE1, EE1_reads)

# Analyze Malagasy annotatations at cluster level
C <- fread("~/dev/figshare-repos/iba/processed_data/MG.v2/cluster_counts.tsv")
C_reads <- row_reads(C)
T2 <- read.delim("~/dev/figshare-repos/iba/processed_data/MG.v2/cluster_taxonomy.tsv")
T2 <- T2[T2$representative==1,]

S2 <- S2[S2$ASV %in% T2$ASV,]
S2_reads <- match_reads(C,C_reads,T2$cluster[match(S2$ASV,T2$ASV)])
annotation_res <- add_results(annotation_res, "sintax", "Madagascar", "clusters", S2, S2_reads)

V2 <- V2[V2$ASV %in% T2$ASV,]
V2_reads <- match_reads(C,C_reads,T2$cluster[match(V2$ASV,T2$ASV)])
annotation_res <- add_results(annotation_res, "vsearch", "Madagascar", "clusters", V2, V2_reads)

E2 <- E2[E2$ASV %in% T2$ASV,]
E2_reads <- match_reads(C,C_reads,T2$cluster[match(E2$ASV,T2$ASV)])
annotation_res <- add_results(annotation_res, "epang", "Madagascar", "clusters", E2, E2_reads)

EE2 <- EE2[EE2$ASV %in% T2$ASV,]
E2_reads <- match_reads(C,C_reads,T2$cluster[match(EE2$ASV,T2$ASV)])
annotation_res <- add_results(annotation_res, "epang_exp", "Madagascar", "clusters", EE2, EE2_reads)

write.table(annotation_res, "../results/annotation_res.tsv", row.names=FALSE, sep="\t")

