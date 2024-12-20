# Echo filter algorithm
echo_filter <- function(counts, matchlist,
                        min_match=84,
                        n_closest=10,
                        min_overlap=0.95,
                        read_ratio_type="max",
                        max_read_ratio=1.0,
                        require_corr=TRUE,
                        max_p_val=0.05) {

  # Reject if less than 2 clusters
  if (nrow(counts)==0)
    return(list(retained_clusters=character(), discarded_clusters=character()))
  if (nrow(counts)==1)
    return(list(retained_clusters=rownames(counts)[1], discarded_clusters=character()))
  
  # Fix colnames if matchlist is "raw"
  colnames(matchlist)[1:3] <- c("asv1","asv2","idty")

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

    # Find sample reads
    reads <- as.numeric(counts[i,])
    n_samples <- sum(reads>0)

    # Find potential parents within distance
    pot_parent_clusters <- matchlist$asv2[(matchlist$asv1==cluster) & (matchlist$idty>=min_match)]

    # Only keep those that have already been processed
    pot_parent_clusters <- pot_parent_clusters[pot_parent_clusters %in% rownames(counts)[1:(i-1)]]

    # Limit number
    if (length(pot_parent_clusters) > n_closest)
      pot_parent_clusters <- pot_parent_clusters[1:n_closest]

    # Cycle over potential parent clusters  
    for (p_clust in pot_parent_clusters) {
        
      p_reads <- as.numeric(counts[p_clust,])
      ix <- reads!=0 & p_reads!=0
      n_samples_overlap <- sum(ix)
      overlap <- n_samples_overlap / n_samples

      if (overlap >= min_overlap) { 
        if (require_corr && n_samples_overlap > 3) {
          fit <- lm(reads[ix]~p_reads[ix])
          k <- as.numeric(fit$coefficients[2])
          p_val <- as.numeric(summary(fit)$coefficients[2,4])
          if (k <= max_read_ratio && k > 0.0 && p_val < max_p_val) {
            discard <- TRUE
            break
          }
        } else {
          if (read_ratio_type=="max") 
            read_ratio <- max(reads[ix] / p_reads[ix]) 
          else 
            read_ratio <- mean(reads[ix] / p_reads[ix]) 
          if (read_ratio <= max_read_ratio) { 
            discard <- TRUE 
            break; 
          } 
        }
      }
    }
    if (discard) 
      discarded_clusters <- c(discarded_clusters, cluster)
    else        
      retained_clusters <- c(retained_clusters, cluster)
  }
  
  #cat ("|\n") # Done processing all clusters

  list(retained_clusters=retained_clusters, discarded_clusters=discarded_clusters)
}

