library(data.table)


# Help function to get prop_samples, mean, and max reads
#   index:  Index of rows to remove (TRUE) or keep (FALSE)
mean_max <- function(dt,index) {

    is_missing <- rowSums(dt[,2:ncol(dt)])==0

    remove_stats <- index & !is_missing
    remove_mean <- rowSums(dt[remove_stats,2:ncol(dt)]) / rowSums(dt[remove_stats,2:ncol(dt)]>0)
    remove_max <- apply(as.matrix(dt[remove_stats,2:ncol(dt)]),1,max)
    remove_prop <- rowMeans(dt[remove_stats,2:ncol(dt)]>0)

    is_missing_keep <- rowSums(dt[!index,2:ncol(dt)])==0
    keep_stats <- !index & !is_missing
    keep_mean <- rowSums(dt[keep_stats,2:ncol(dt)]) / rowSums(dt[keep_stats,2:ncol(dt)]>0)
    keep_max <- apply(as.matrix(dt[keep_stats,2:ncol(dt)]),1,max)
    keep_prop <- rowMeans(dt[keep_stats,2:ncol(dt)]>0)

    list(remove_mean=remove_mean,
         remove_max=remove_max,
         remove_prop=remove_prop,
         remove_missing=sum(is_missing & index),
         keep_mean=keep_mean,
         keep_max=keep_max,
         keep_prop=keep_prop,
         keep_missing=sum(is_missing & !index))
}


# Function to identify control clusters
#
#   counts:     cluster counts for all samples (including controls)
#   samples:    names of samples
#   controls:   names of controls
#   cutoff:     threshold for removing control clusters
identify_control_clusters <- function(counts, taxonomy, samples, controls, cutoff=0.05) {

    # Identify clusters that appear in controls above cutoff

    # Extract sample and control reads
    idx <- which(colnames(counts) %in% samples)
    idx <- c(1,idx) # keep cluster name
    sample_counts <- counts[,..idx]

    idx <- which(colnames(counts) %in% controls)
    idx <- c(1,idx) # keep cluster name
    control_counts <- counts[,..idx]

    include_rows <- rowSums(control_counts[, 2:ncol(control_counts)])>0
    control_counts <- control_counts[include_rows,]
    sample_counts <- sample_counts[include_rows,]

    # Output some diagnostics
    tot_clusters <- nrow(control_counts)
    prop_controls <- rowMeans(control_counts[,2:ncol(control_counts)]>0)

    cat("There are", tot_clusters, "clusters in the control samples\n")
    cat("Of these,", sum(prop_controls>cutoff), "occur in more than", cutoff, "of samples\n")

    remove_clusters <- control_counts$cluster[prop_controls>cutoff]
    
    res <- mean_max(control_counts, prop_controls>cutoff)
    cat("Reads in controls of removed clusters:\n")
    remove_tax <- data.frame(list(cluster=remove_clusters,
                                  prop_controls=res$remove_prop,
                                  mean_reads=res$remove_mean,
                                  max_reads=res$remove_max,
                                  Genus=taxonomy$Genus[match(remove_clusters,taxonomy$cluster)],
                                  Species=taxonomy$Species[match(remove_clusters,taxonomy$cluster)],
                                  BOLD_bin=taxonomy$BOLD_bin[match(remove_clusters,taxonomy$cluster)]))
    print(remove_tax)
    cat("Summary of reads in controls of kept clusters:\n")
    cat("prop_samples:\n")
    print(summary(res$keep_prop))
    cat("mean:\n")
    print(summary(res$keep_mean))
    cat("max:\n")
    print(summary(res$keep_max))

    prop_samples <- rowMeans(sample_counts[,2:ncol(sample_counts)]>0)
    res <- mean_max(sample_counts, prop_controls>cutoff)
    cat("Reads in samples of removed clusters:\n")
    remove_tax <- data.frame(list(cluster=remove_clusters,
                                  prop_samples=res$remove_prop,
                                  mean_reads=res$remove_mean,
                                  max_reads=res$remove_max,
                                  Genus=taxonomy$Genus[match(remove_clusters,taxonomy$cluster)],
                                  Species=taxonomy$Species[match(remove_clusters,taxonomy$cluster)],
                                  BOLD_bin=taxonomy$BOLD_bin[match(remove_clusters,taxonomy$cluster)]))
    print(remove_tax)
    cat("Summary of reads in samples of kept clusters:\n")
    cat("prop_samples:\n")
    print(summary(res$keep_prop))
    cat("mean:\n")
    print(summary(res$keep_mean))
    cat("max:\n")
    print(summary(res$keep_max))

    remove_clusters
}


# Function to remove control_clusters
#
#   filtered_counts:    filtered cluster counts for samples, col 1 should be cluster name
#   all_cluster_counts: cluster counts for all samples
#   samples:            names of samples
#   controls:           names of controls
#   cutoff:             threshold for removing control clusters
remove_control_clusters <- function(filtered_counts, all_cluster_counts, taxonomy, samples, controls, cutoff=0.05) {

    remove_clusters <- identify_control_clusters(all_cluster_counts, taxonomy, samples, controls, cutoff)

    filtered_counts[!filtered_counts$cluster %in% remove_clusters,]
}


# Function for identifying spikeins
identify_spikes <- function(counts, spikein_samples, taxonomy, cutoff=0.8) {

    # Get a rep asv taxonomy in case a complete cluster taxonomy is provided
    if ("representative" %in% colnames(taxonomy))
        taxonomy <- taxonomy[taxonomy$representative==1,]

    # Get counts for the samples containing spikeins
    # Account for counts being either data.frame or data.table
    idx <- which(colnames(counts) %in% spikein_samples)
    idx <- c(1,idx)
    if (class(counts)[1]=="data.table")
        counts <- counts[,..idx]
     else
        counts <- counts[,idx]

    # Identify spikeins
    prop_samples <- rowMeans(counts[,2:ncol(counts)]>0)
    spikein_candidates <- counts$cluster[prop_samples>cutoff]
    spikein_clusters <- spikein_candidates[taxonomy$Class[match(spikein_candidates,taxonomy$cluster)]=="Insecta"]
    cat("Found", length(spikein_clusters), "spikein clusters:\n")
    spike_tax <- taxonomy[taxonomy$cluster %in% spikein_clusters,c("cluster","Genus","Species","BOLD_bin")]
    spike_tax$prop_samples <- rowMeans(counts[match(spike_tax$cluster,counts$cluster),2:ncol(counts)]>0)
    #print(spike_tax)

    # Return clusters
    spikein_clusters
}


# Function for removing spikes
remove_spikes <- function(counts, spikein_samples, taxonomy, cutoff=0.8) {

    spikein_clusters <- identify_spikes(counts, spikein_samples, taxonomy)

    counts[!(counts$cluster %in% spikein_clusters),]
}

