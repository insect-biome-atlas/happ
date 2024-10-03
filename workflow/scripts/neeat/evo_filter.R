# Evo filter algorithm
evo_filter <- function(counts, evodistlist,
                        min_match=84,
                        n_closest=10,
                        require_overlap=TRUE,
                        min_overlap=0.95,
                        dist_type="wdadn",
                        dist_threshold=1.0,
                        remove_exotics=FALSE) {
  
  # Reject if less than 2 clusters
  if (nrow(counts)==0)
    return(list(retained_clusters=character(), discarded_clusters=character()))
  if (nrow(counts)==1)
    return(list(retained_clusters=rownames(counts)[1], discarded_clusters=character()))

  # Get total number of reads of each cluster
  cat("Getting total number of samples and reads for each cluster\n")
  tot_samples <- rowSums(counts>0)
  tot_reads <- rowSums(counts)

  # Order clusters 
  cat ("Ordering clusters\n")
  ix <- order(tot_samples, tot_reads, decreasing=TRUE)
  counts <- counts[ix,]
  tot_samples <- tot_samples[ix]
  tot_reads <- tot_reads[ix]

  # Cycle over clusters

  cat (paste0("Processing ", nrow(counts), " clusters...\n"))
  #cat ("0%                                                100%\n")
  #cat ("|")
  print_interval <- ceiling(nrow(counts) / 50)

  retained_clusters <- rownames(counts)[1]
  discarded_clusters <- character()
  pb <- txtProgressBar(min = 0, max = nrow(counts)-1, style = 3, width = 50)
  for (i in 2:nrow(counts)) {

    # Print progress
    #if (i %% print_interval == 0)
    #  cat("-")
    setTxtProgressBar(pb, i)
    # Find cluster name
    cluster <- rownames(counts)[i]
    discard <- FALSE

    # Find sample reads and n_samples
    reads <- as.numeric(counts[i,])
    n_samples <- sum(reads>0)

    # Find potential parents within distance
    pot_parent_clusters <- evodistlist$asv2[(evodistlist$asv1==cluster) & (evodistlist$idty>=min_match)]
    if (remove_exotics && length(pot_parent_clusters)==0)
      exotic <- TRUE
    else
      exotic <- FALSE

    # Only keep those that have passed the filter
    pot_parent_clusters <- pot_parent_clusters[pot_parent_clusters %in% retained_clusters]

    # Limit number
    if (length(pot_parent_clusters) > n_closest)
      pot_parent_clusters <- pot_parent_clusters[1:n_closest]

    # Find match values
    local_evodistlist <- evodistlist[(evodistlist$asv1==cluster) & (evodistlist$asv2 %in% pot_parent_clusters),]

    # Cycle over potential parent clusters  
    for (p_clust in pot_parent_clusters) {

      if (require_overlap) {
        p_reads <- as.numeric(counts[p_clust,])
        n_samples_overlap <- sum(reads!=0 & p_reads!=0)
        overlap_frac <- n_samples_overlap / n_samples
        if (overlap_frac < min_overlap)
          next
      }

      p_dist <- local_evodistlist$pdist[match(p_clust,local_evodistlist$asv2)]
      if (dist_type=="dadn") {
        e_dist <- local_evodistlist$dadn[match(p_clust,local_evodistlist$asv2)]
      } else {
        e_dist <- local_evodistlist$wdadn[match(p_clust,local_evodistlist$asv2)]
      }
      # p_dist > 0 by definition (ASVs are uniqe)
      # we need not worry about division with zero here
      if (e_dist/p_dist > dist_threshold) {
        discard <- TRUE
        break
      }
    }
    if (discard || exotic) 
      discarded_clusters <- c(discarded_clusters, cluster)
    else        
      retained_clusters <- c(retained_clusters, cluster)
  }
  
  #cat ("|\n") # Done processing all clusters

  list(retained_clusters=retained_clusters, discarded_clusters=discarded_clusters)
}

