sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)
source(snakemake@params$functions)
library(data.table)

## INPUT
counts_file <- snakemake@input$counts
taxonomy_file <- snakemake@input$taxonomy

## PARAMS
meta_file <- snakemake@params$meta
sample_id_col <- snakemake@params$sample_id_col
sample_type_col <- snakemake@params$sample_type_col
sample_val <- snakemake@params$sample_val

# Read in data
cat("Processing data\n")
#meta <- get_se_meta()

if ((meta_file!=FALSE) && file.exists(meta_file)) {
    meta <- read.delim(meta_file)
    rownames(meta) <- meta[, 1]
    if (sample_id_col != "") {
        rownames(meta) <- meta[, sample_id_col]
    }
    if ((sample_type_col!="") && (length(sample_val)>0)) {
        meta <- meta[meta[,sample_type_col]%in%sample_val,]
    }
    samples <- rownames(meta)
} else {
    meta <- FALSE
    samples <- c()
}

counts <- fread(counts_file, sep="\t")
counts <- as.data.frame(counts)
taxonomy <- read.delim(taxonomy_file)

cat("Generating raw counts table\n")
# extract counts for samples
if (is.null(samples)) {
    index <- 2:ncol(counts)
} else {
    index <- match(samples,colnames(counts))
}
index <- c(1, index[!is.na(index)])
counts <- counts[,index]
# extract counts with rowSum > 0
include_rows <- rowSums(counts[,2:ncol(counts)])>0
counts <- counts[include_rows,]
# extract counts with colSum > 0
include_cols <- c(TRUE,as.logical(colSums(counts[,2:ncol(counts)])>0))
counts <- counts[,include_cols]
#write.table(counts,"../bigdata/cluster_counts.tsv", row.names=FALSE, sep="\t")
write.table(counts, snakemake@output$cluster_counts, row.names=FALSE, sep="\t")


cat("Generating calibrated counts table\n")

# identify spikeins
# Require that the metadata file has a column "spikein_sample" that is TRUE for spike-in samples
#spikein_samples <- meta$sampleID_NGI[meta$dataset %in% c("CO1_lysate_2019_SE","CO1_homogenate_2019_SE")]
if ((!is.null(samples)) && ("spikein_sample" %in% colnames(meta))) {
    meta$spikein_sample <- as.logical(meta$spikein_sample)
    spikein_samples <- rownames(meta[meta$spikein_sample==TRUE,])
    spikein_clusters <- identify_spikes(counts, spikein_samples, taxonomy[taxonomy$representative==1,])
} else {
    spikein_clusters <- c()
}

# calibrate
# Note that there are samples without spike-ins; we do not recalibrate these.
# If no meta file is supplied, the calibrated counts will be the same as the raw counts.
ccounts <- counts
if (length(spikein_clusters) > 0) {
    cat("Calibrating counts\n")
    spike_counts <- colSums(ccounts[ccounts$cluster %in% spikein_clusters,2:ncol(ccounts)])
    correction <- spike_counts / mean(spike_counts[spike_counts!=0])
    correction[spike_counts==0] <- 1.0
    for (i in 2:ncol(ccounts)) {
        if (i%%100==0)
            cat("Processing col ", i, " (", round(100*i/ncol(counts)), "%)\n",sep="")
        ccounts[,i] <- ceiling(ccounts[,i] / correction[i-1])
    }
} else {
    cat("No spike-ins found, not calibrating counts\n")
}
#write.table(ccounts, "../bigdata/calibrated_cluster_counts.tsv", sep="\t",row.names=FALSE)
write.table(ccounts, snakemake@output$calibrated_counts, sep="\t",row.names=FALSE)
rm(ccounts)


cat("Generating proportional (of total) counts table\n")

pcounts <- counts
sample_counts <- colSums(pcounts[,2:ncol(pcounts)])
for (i in 2:ncol(pcounts)) {
    if (i%%100==0)
        cat("Processing col ", i, " (", round(100*i/ncol(pcounts)), "%)\n",sep="")
    if (sample_counts[i-1] != 0)
        pcounts[,i] <- pcounts[,i] / sample_counts[i-1]
}
#write.table(pcounts,"../bigdata/tot_proportional_cluster_counts.tsv",sep="\t",row.names=FALSE)
write.table(pcounts, snakemake@output$tot_proportional_counts, sep="\t",row.names=FALSE)

if (length(spikein_clusters) > 0) {
    pcounts <- counts
    cat("Generating proportional (of sample reads without spikeins) counts table\n")
    for (i in 2:ncol(pcounts)) {
        if (i%%100==0)
            cat("Processing col ", i, " (", round(100*i/ncol(pcounts)), "%)\n",sep="")
        if (sample_counts[i-1] != 0)
            pcounts[,i] <- pcounts[,i] / (sample_counts[i-1] - spike_counts[i-1])
    }
} 
#write.table(pcounts,"../bigdata/sample_proportional_cluster_counts.tsv",sep="\t",row.names=FALSE)
write.table(pcounts, snakemake@output$sample_proportional_counts, sep="\t",row.names=FALSE)

sink()