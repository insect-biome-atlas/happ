source("~/dev/ms-repos-iba/utils/spikes_controls_fxns.R")
source("get_data_fxns.R")

# Read in Sweden data
cat("Processing Sweden data\n")
meta <- get_se_meta()
counts <- get_se_cluster_counts_df()
taxonomy <- get_se_cluster_taxonomy()

samples <- meta$sampleID_NGI[meta$lab_sample_type=="sample" & meta$sequencing_status=="sequencing successful"]

cat("Generating raw counts table\n")
index <- match(samples,colnames(counts))
index <- c(1, index[!is.na(index)])
counts <- counts[,index]
include_rows <- rowSums(counts[,2:ncol(counts)])>0
counts <- counts[include_rows,]
include_cols <- c(TRUE,as.logical(colSums(counts[,2:ncol(counts)])>0))
counts <- counts[,include_cols]
write.table(counts,"../bigdata/cluster_counts.tsv", row.names=FALSE, sep="\t")


cat("Generating calibrated counts table\n")

# identify spikeins
spikein_samples <- meta$sampleID_NGI[meta$dataset %in% c("CO1_lysate_2019_SE","CO1_homogenate_2019_SE")]
spikein_clusters <- identify_spikes(counts, spikein_samples, taxonomy[taxonomy$representative==1,])

# calibrate
# Note that there are samples without spike-ins; we do not recalibrate these.
ccounts <- counts
if (length(spikein_clusters) > 0) {
    spike_counts <- colSums(ccounts[ccounts$cluster %in% spikein_clusters,2:ncol(ccounts)])
    correction <- spike_counts / mean(spike_counts[spike_counts!=0])
    correction[spike_counts==0] <- 1.0
    for (i in 2:ncol(ccounts)) {
        if (i%%10==0)
            cat("Processing col ", i, " (", round(100*i/ncol(counts)), "%)\n",sep="")
        ccounts[,i] <- ceiling(ccounts[,i] / correction[i-1])
    }
}
write.table(ccounts, "../bigdata/calibrated_cluster_counts.tsv", sep="\t",row.names=FALSE)
rm(ccounts)


cat("Generating proportional (of total) counts table\n")

pcounts <- counts
sample_counts <- colSums(pcounts[,2:ncol(pcounts)])
for (i in 2:ncol(pcounts)) {
    if (i%%10==0)
        cat("Processing col ", i, " (", round(100*i/ncol(pcounts)), "%)\n",sep="")
    if (sample_counts[i-1] != 0)
        pcounts[,i] <- pcounts[,i] / sample_counts[i-1]
}
write.table(pcounts,"../bigdata/tot_proportional_cluster_counts.tsv",sep="\t",row.names=FALSE)


cat("Generating proportional (of sample reads without spikeins) counts table\n")

pcounts <- counts
for (i in 2:ncol(pcounts)) {
    if (i%%10==0)
        cat("Processing col ", i, " (", round(100*i/ncol(pcounts)), "%)\n",sep="")
    if (sample_counts[i-1] != 0)
        pcounts[,i] <- pcounts[,i] / (sample_counts[i-1] - spike_counts[i-1])
}
write.table(pcounts,"../bigdata/sample_proportional_cluster_counts.tsv",sep="\t",row.names=FALSE)

