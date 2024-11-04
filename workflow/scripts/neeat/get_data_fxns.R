# Functions for getting original data. Modify the paths to fit your system.
get_se_meta <- function() { read.delim("~/dev/figshare-repos/iba/raw_data/CO1_sequencing_metadata_SE.tsv") }

get_se_cluster_taxonomy <- function() { read.delim("~/dev/figshare-repos/iba/processed_data/SE.v2/cluster_taxonomy.tsv") }
get_se_cluster_rep_taxonomy <- function() { T<- get_se_cluster_taxonomy(); T[T$representative==1,] }

#get_se_cluster_counts_df <- function() { read.delim("~/dev/figshare-repos/iba/processed_data/SE.v2/cluster_counts.tsv") }
get_se_cluster_counts_dt <- function() { require(data.table); fread("~/dev/figshare-repos/iba/processed_data/SE.v2/cluster_counts.tsv") }

get_se_cluster_taxonomy_samples_hexapoda <- function() {

    meta <- get_se_meta()

    samples <- meta$sampleID_NGI[meta$lab_sample_type=="sample" & meta$sequencing_status=="sequencing successful"]

    counts <- get_se_cluster_counts_dt()

    index <- match(samples,colnames(counts))
    index <- c(1, index[!is.na(index)])
    counts <- counts[,..index]
    include_rows <- rowSums(counts[,2:ncol(counts)])>0
    counts <- counts[include_rows,]
    include_cols <- c(TRUE,as.logical(colSums(counts[,2:ncol(counts)])>0))
    counts <- counts[,..include_cols]

    taxonomy <- get_se_cluster_taxonomy()

    hex_classes <- c("Insecta","Diplura","Protura","Collembola")
    taxonomy <- taxonomy[taxonomy$Class %in% hex_classes,]

    taxonomy[taxonomy$cluster %in% counts$cluster,]
}

