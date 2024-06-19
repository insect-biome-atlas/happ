evo_dist <- function(x, y, start_pos=1) {

    x <- toupper(x)
    y <- toupper(y)
    
    dist <- 0

    for (i in seq(from=start_pos, to=nchar(x)-2, by=3)) {

        codon1 <- substring(x,i,i+2)
        codon2 <- substring(y,i,i+2)
        
        if (!grepl("-", codon1) && !grepl("-", codon2)) {
            dist <- dist + codon_dist[codon1,codon2]
        } else if ((grepl("---", codon1) && !grepl("-", codon2)) ||
                   (grepl("---", codon2) && !grepl("-", codon1))) {
          dist <- dist + 100
        }
    }

    dist / (nchar(x)-start_pos+1)
}

dadn_dist <- function(x, y, start_pos=1) {
  
  x <- toupper(x)
  y <- toupper(y)
  
  da <- 0
  dn <- 0
  
  for (i in seq(from=start_pos, to=nchar(x)-2, by=3)) {
    
    codon1 <- substring(x,i,i+2)
    codon2 <- substring(y,i,i+2)

    if (!grepl("-", codon1) && !grepl("-", codon2)) {
      da <- da + codon_aa_dist[codon1,codon2]
      dn <- dn + codon_nuc_dist[codon1,codon2]
    } else if ((grepl("---", codon1) && !grepl("-", codon2)) ||
               (grepl("---", codon2) && !grepl("-", codon1))) {
      da <- da + 1
    }
      
  }
  
  da / dn
}


wdadn_dist <- function(x, y, start_pos=1) {

  x <- toupper(x)
  y <- toupper(y)
  
  da <- 0
  dn <- 0
  
  for (i in seq(from=start_pos, to=nchar(x)-2, by=3)) {
    
    codon1 <- substring(x,i,i+2)
    codon2 <- substring(y,i,i+2)
    
    if (!grepl("-", codon1) && !grepl("-", codon2)) {
      da <- da + codon_waa_dist[codon1,codon2]
      dn <- dn + codon_nuc_dist[codon1,codon2]
    } else if ((grepl("---", codon1) && !grepl("-", codon2)) ||
               (grepl("---", codon2) && !grepl("-", codon1))) {
      da <- da + 1
    }
  }
  
  da / dn
}

abundance_filter <- function(order_name, taxonomy, counts, output, save_to_file=TRUE, threshold=3) {

  # Subsetting taxonomy data
  cat("Subsetting taxonomy data\n")
  taxonomy <- taxonomy[taxonomy$Order==order_name & taxonomy$representative,]
  
  # Get total number of reads of each cluster
  cat("Getting total number of reads for each cluster\n")
  rowsums <- data.frame(rowSums(counts[, 2:ncol(counts)]))
  rowsums$cluster <- counts$cluster
  colnames(rowsums) <- c("n_reads", "cluster")

  # Order clusters in terms of total number of reads
  cat ("Ordering clusters in terms of reads\n")
  ix <- order(rowsums$n_reads, decreasing=TRUE)
  clust_sum <- rowsums[ix,]

  # Get results
  res <- data.frame(list(cluster=clust_sum$cluster, numt=clust_sum$n_reads < threshold))
  
  if (save_to_file) {
    write.table(res,output, sep="\t", row.names=FALSE, quote=FALSE)
  }
  res
}

combined_filter_neighbors <- function(fasta, taxonomy, counts, output, spikeins=NULL,
                        distance="wdadn", threshold=1.2, relative=TRUE,
                        max_corr_coeff=0.010, max_p_val=0.05, max_read_ratio=0.10,
                        max_singletons=1, max_singleton_reads=20, save_to_file=TRUE, n_closest=5) {
  # Read in fasta alignment 
  cat("Reading alignment\n")
  seqs <- read.alignment(fasta,format="fasta", whole.header = FALSE)
  
  # Name the sequences after the cluster instead of the ASV
  sequence_names <- taxonomy$cluster[match(seqs$nam,taxonomy$ASV)]
  
  
  # Compute distance matrix and label rows and cols in matrix with cluster name
  cat("Computing distance matrix\n")
  dist_matr <- as.matrix(dist.alignment(seqs))
  rownames(dist_matr) <- taxonomy$cluster[match(rownames(dist_matr),taxonomy$ASV)]
  colnames(dist_matr) <- taxonomy$cluster[match(colnames(dist_matr),taxonomy$ASV)]
  
  # Get and process counts data
  
  if (length(spikeins) > 0) {
    cat ("Extracting and calibrating counts data\n")
    spikein_reads <- colSums(counts[counts$cluster %in% spikeins,2:ncol(counts)])
    ix <- spikein_reads != 0
    spikein_reads <- spikein_reads[ix]
    mean_spikein_reads <- mean(spikein_reads) 
    counts <- counts[,c(1,(2:ncol(counts))[ix])]

    # Extract relevant rows in counts data
    counts <- counts[counts$cluster %in% rownames(dist_matr),]
    
    # Calibrate data using spikeins
    # Note the difference in index in spikein_reads (first element of spikein_reads is column 2 in counts)
    for (i in 2:ncol(counts))
      counts[,i] <- ceiling(mean_spikein_reads * counts[,i]/spikein_reads[i-1])
  }
  
  # Get total number of reads of each cluster
  cat("Getting total number of reads for each cluster\n")
  # using rowSum here since it's faster
  tot_reads <- as.vector(rowSums(counts[, 2:ncol(counts)]))
  
  # Order clusters in terms of total number of reads
  cat ("Ordering clusters in terms of reads\n")
  ix <- order(tot_reads, decreasing=TRUE)
  counts <- counts[ix,]
  tot_reads <- tot_reads[ix]
  
  # Prepare data table for results 
  reads <- as.numeric(counts[1,2:ncol(counts)])
  num_rows <- length(unique(counts$cluster))-1
  num_columns <- 5

  res <- data.table(matrix(NA, nrow = num_rows, ncol = num_columns))
  colnames(res) <- c(
                  "cluster",
                  "n_samples",
                  "n_reads",
                  "numt",
                  "reason")

    first_cluster <- list(counts$cluster[1],
                  sum(reads != 0),
                  tot_reads[1],
                  FALSE,
                   "")

  res <- rbindlist(list(res, first_cluster))
  row_index <- 1

  # Cycle over clusters (the most abundant cluster is not a numt)
  cat ("Analyzing clusters\n")
  co1_clusters <- counts$cluster[1]
  numt_clusters <- character()
  n_iter <- length(2:nrow(counts))
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = n_iter, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar
  for (i in 2:nrow(counts)) {    
    # Find cluster name
    cluster <- counts$cluster[i]
    # Find sequence
    s1 <- seqs$seq[[which(sequence_names==cluster)]]

    # Test whether it is a numt
    isNumt <- FALSE
    reason <- ""
    
    # Find sample reads
    reads <- as.numeric(counts[counts$cluster==cluster,2:ncol(counts)])
    
    # Find samples where cluster exist: cluster_samples
    cluster_samples_indeces <- which(counts[counts$cluster == cluster, -1] > 0)
    cluster_samples <- colnames(counts)[cluster_samples_indeces+1]

  ###################################
    if (sum(reads != 0) > 3) { # Branch 1: If cluster in > 3 samples
      ## Branch 1.1 "p_deviation"

      # Subset counts table on co1_clusters and cluster_samples
      co1clusters_row_indices <- which(counts$cluster %in% co1_clusters)
      subset_counts <- counts[co1clusters_row_indices, c("cluster", cluster_samples)]

      # Choose/Filter the CO1 clusters that exist in all but max 1 of cluster_samples
      cluster_presence <- rowMeans(subset_counts[, -1] > 0)
      threshold_presence <- (length(cluster_samples)-1)/length(cluster_samples) 
      filtered_co1clusters <- subset_counts$cluster[cluster_presence >= threshold_presence] 

      # Rank these CO1-clusters by distance to cluster    
      subset_matrix <- as.matrix(dist_matr[cluster,filtered_co1clusters]) #subset to found CO1_clusters and cluster
      colnames(subset_matrix) <- cluster
      column_values <- (as.numeric(subset_matrix[,1]))^2 #p-dist
      sorted_indices <- order(column_values)
      most_closely_related_co1clusters <- rownames(subset_matrix)[sorted_indices]  
      
      # Choose n_closest most closely related filtered co1 clusters
      if (length(filtered_co1clusters) < n_closest) {
        pot_parent_clusters <- filtered_co1clusters
      } else 
      {
        pot_parent_clusters <- most_closely_related_co1clusters[1:n_closest]
      }

      # Go through these chosen CO1-clusters, "pot_parent_clusters", calculate dist_value
      if ((length(pot_parent_clusters)) > 0) {
        for (p_clust in pot_parent_clusters) {
          d <- dist_matr[cluster, p_clust]
          p_dist <- d^2
          s2 <- seqs$seq[[which(sequence_names==p_clust)]]
          e_dist <- evo_dist(s1,s2)
          dadn <- dadn_dist(s1,s2)
          wdadn <- wdadn_dist(s1,s2)  

          # Find the critical value to filter on
          if (distance=="wdadn") dist_value <- wdadn
          else if (distance=="dadn") dist_value <- dadn
          else dist_value <- e_dist
          if (relative)
            dist_value <- dist_value / p_dist

          if (dist_value > threshold) { # ID as numt if strange sequence at aa-level
            isNumt <- TRUE
            reason <- "p_deviation"
          }
        }
      }

      if (!isNumt) {
      ## Branch 1.2 "p_correlation"
        if ((length(pot_parent_clusters)) > 0) {
          for (p_clust in pot_parent_clusters) {
            p_reads <- as.numeric(counts[counts$cluster==p_clust,2:ncol(counts)])
            ix <- reads!=0 & p_reads!=0
            n_samples_overlap <- sum(ix)
            if (n_samples_overlap > 3) {
              fit <- lm(reads[ix]~p_reads[ix])
              corr_coeff <- as.numeric(fit$coefficients[2])
              p_val <- as.numeric(summary(fit)$coefficients[2,4])
              if (corr_coeff < max_corr_coeff && corr_coeff > 0.0 && p_val < max_p_val) {
                isNumt <- TRUE
                reason <- "p_correlation"
              }
            }
          }
        }
      }
    }

    ###
    else { # Branch 2: If cluster in < 3 samples

    # Branch 2.1 "small_sample"

    # Choose all the CO1 clusters 
    filtered_co1clusters <- co1_clusters

    # Rank these CO1-clusters by distance to cluster    
    subset_matrix <- as.matrix(dist_matr[cluster, filtered_co1clusters]) #subset to found CO1_clusters and cluster
    colnames(subset_matrix) <- cluster
    column_values <- (as.numeric(subset_matrix[,1]))^2 #p-dist
    sorted_indices <- order(column_values)
    most_closely_related_co1clusters <- rownames(subset_matrix)[sorted_indices]   

    # Choose n_closest*2 most closely related filtered co1 clusters
    if (length(filtered_co1clusters) < (n_closest*2)) {
      pot_parent_clusters <- filtered_co1clusters
    }
    else {
      pot_parent_clusters <- most_closely_related_co1clusters[1:(n_closest*2)]
    }

    # Go through these chosen CO1-clusters, "pot_parent_clusters"
    if ((length(pot_parent_clusters)) > 0) {
      for (p_clust in pot_parent_clusters) {
        d <- dist_matr[cluster,p_clust]
        p_dist <- d^2
        s2 <- seqs$seq[[which(sequence_names==p_clust)]]
        e_dist <- evo_dist(s1,s2)
        dadn <- dadn_dist(s1,s2)
        wdadn <- wdadn_dist(s1, s2) 

        # Find the critical value to filter on
        if (distance=="wdadn") dist_value <- wdadn
        else if (distance=="dadn") dist_value <- dadn
        else dist_value <- e_dist
        if (relative)
          dist_value <- dist_value / p_dist
    
        p_reads <- as.numeric(counts[counts$cluster==p_clust,2:ncol(counts)])
        ix <- reads!=0 & p_reads!=0
        n_samples_overlap <- sum(ix)

        if (n_samples_overlap == 0) {
          ratio_max <- 0
        }
        else {
          ratio_max <- max(reads[ix] / p_reads[ix])
        }
      
        singletons <- sum(reads!=0) - sum(ix)
        if (singletons > 0)
          singleton_reads <- median(reads[reads!=0 & p_reads==0])
        else
          singleton_reads <- NA

        if (ratio_max < max_read_ratio && dist_value > threshold && (singletons==0 || (singletons<=max_singletons && singleton_reads<=max_singleton_reads))) {
          isNumt <- TRUE
          reason <- "small_sample"
        }
      }
    }
  }

  ###################################
  ## Results
  if (!isNumt) 
    co1_clusters <- c(co1_clusters, cluster)
  else        
    numt_clusters <- c(numt_clusters, cluster)

  new_row <- data.table(cluster, sum(reads!=0), sum(reads),
                        isNumt, reason)

  res[row_index, colnames(res) := new_row]
  row_index <- row_index + 1
  setTxtProgressBar(pb, i)
  }

  close(pb)

  res <- res[order(res$n_reads, decreasing = TRUE),]

  #####
  if (save_to_file) {
    write.table(res,output, sep="\t", row.names=FALSE, quote=FALSE)
  }
  ###

  #Return results
  res

}

evaluate_res <- function(order_name, order_clusters, taxonomy, numt_data, trusted_file="", outfile) {
  T <- taxonomy
  cat ("Reading taxonomy and cluster data\n")
  numt_clusters <- numt_data$cluster[numt_data$numt]
  if (trusted_file != "") {
    cat ("Reading trusted data\n")
    trusted_asvs <- read.table(trusted_file, header=FALSE)
    colnames(trusted_asvs) <- c("ASV")
    trusted_clusters <- unique(T$cluster[match(trusted_asvs$ASV,T$ASV)])
    trusted <- sum(order_clusters$cluster %in% trusted_clusters)
    removed_trusted <- sum(numt_clusters %in% trusted_clusters)
    cat ("False positive rate (trusted clusters removed): ", removed_trusted, "of", trusted, "=", removed_trusted/trusted, "\n", file=outfile, append=TRUE)
  }

  cat (order_name, "has", nrow(order_clusters), "clusters in total\n", file=outfile, append=TRUE)
  cat (sum(grepl(" ",order_clusters$Species)), "clusters have rep asv annotated to species\n", file=outfile, append=TRUE)
  cat ("There are", length(unique(order_clusters$Species[grepl(" ",order_clusters$Species)])), "unique species annotations\n", file=outfile, append=TRUE)
  cat (sum(grepl("BOL",order_clusters$BOLD_bin)), "clusters have rep asv annotated to BOLD bin\n", file=outfile, append=TRUE)
  cat ("There are", length(unique(order_clusters$BOLD_bin[grepl("BOL",order_clusters$BOLD_bin)])), "unique BOLD bin annotations\n", file=outfile, append=TRUE)
  cat (length(numt_clusters), "numts removed\n", file=outfile, append=TRUE)

  unique_species <- length(unique(order_clusters$Species[grepl(" ", order_clusters$Species)]))
  remaining_unique_species <- length(unique(order_clusters$Species[grepl(" ", order_clusters$Species) & !(order_clusters$cluster %in% numt_clusters)]))
  cat ("False positive rate (unique species removed): ", unique_species-remaining_unique_species, "of", unique_species, "=", (unique_species-remaining_unique_species)/unique_species, "\n", file=outfile, append=TRUE)
  
  unique_bins <- length(unique(order_clusters$BOLD_bin[grepl("BOL", order_clusters$BOLD_bin)]))
  remaining_unique_bins <- length(unique(order_clusters$BOLD_bin[grepl("BOL", order_clusters$BOLD_bin) & !(order_clusters$cluster %in% numt_clusters)]))
  cat ("False positive rate (unique BOLD bins removed): ", unique_bins-remaining_unique_bins, "of", unique_bins, "=", (unique_bins-remaining_unique_bins)/unique_bins, "\n", file=outfile, append=TRUE)
  
  dups <- sum(duplicated(order_clusters$Species[grepl(" ", order_clusters$Species)]))
  remaining_dups <- sum(duplicated(order_clusters$Species[grepl(" ", order_clusters$Species) & !(order_clusters$cluster %in% numt_clusters)]))
  cat ("False negative rate (duplicated species records remaining): ", remaining_dups, "of", dups, "=", remaining_dups/dups, "\n", file=outfile, append=TRUE)

  dups <- sum(duplicated(order_clusters$BOLD_bin[grepl("BOL", order_clusters$BOLD_bin)]))
  BOL_bins <- order_clusters[grepl("BOL", order_clusters$BOLD_bin), ]
  non_numt_BOL_bins <- BOL_bins[!BOL_bins$cluster %in% numt_clusters, ]
  remaining_dups <- sum(duplicated(non_numt_BOL_bins$BOLD_bin))
  cat ("False negative rate (duplicated BOLD bin records remaining): ", remaining_dups, "of", dups, "=", remaining_dups/dups, "\n", file=outfile, append=TRUE)

  dups <- length(unique(order_clusters$Species[grepl(" ", order_clusters$Species) & duplicated(order_clusters$Species)]))
  correct_species <- order_clusters[grepl(" ", order_clusters$Species), ]
  non_numt_correct_species <- correct_species[!correct_species$cluster %in% numt_clusters, ]
  remaining_dups <- sum(table(non_numt_correct_species$Species) > 1)
  cat ("False negative rate (duplicated species remaining): ", remaining_dups, "of", dups, "=", remaining_dups/dups, "\n", file=outfile, append=TRUE)
  
  dups <- length(unique(order_clusters$BOLD_bin[grepl("BOL", order_clusters$BOLD_bin) & duplicated(order_clusters$BOLD_bin)]))
  BOL_bins <- order_clusters[grepl("BOL", order_clusters$BOLD_bin), ]
  non_numt_BOL_bins <- BOL_bins[!BOL_bins$cluster %in% numt_clusters, ]
  remaining_dups <- length(unique(non_numt_BOL_bins$BOLD_bin[grepl("BOL", non_numt_BOL_bins$BOLD_bin) & duplicated(non_numt_BOL_bins$BOLD_bin)]))
  cat ("False negative rate (duplicated BOLD bins remaining): ", remaining_dups, "of", dups, "=", remaining_dups/dups, "\n", file=outfile, append=TRUE)
}

evaluate_res_details <- function(order_name, order_clusters, taxonomy, numt_data, trusted_file, outfile) {
  
  evaluate_res(order_name, order_clusters, taxonomy, numt_data, trusted_file, outfile)
  
  numt_clusters <- numt_data$cluster[numt_data$numt]

  if ("Species_resolved" %in% colnames(order_clusters)) {
    unique_species <- length(unique(order_clusters$Species_resolved[grepl(" ", order_clusters$Species_resolved)]))
    remaining_unique_species <- length(unique(order_clusters$Species_resolved[grepl(" ", order_clusters$Species_resolved) & !(order_clusters$cluster %in% numt_clusters)]))
    cat ("False positive rate (unique resolved species removed): ", unique_species-remaining_unique_species, "of", unique_species, "=", (unique_species-remaining_unique_species)/unique_species, "\n", file=outfile, append=TRUE)

    dups <- sum(duplicated(order_clusters$Species_resolved[grepl(" ", order_clusters$Species_resolved)]))
    remaining_dups <- sum(duplicated(order_clusters$Species_resolved[grepl(" ", order_clusters$Species_resolved) & !(order_clusters$cluster %in% numt_clusters)]))
    cat ("False negative rate (duplicated resolved species records remaining): ", remaining_dups, "of", dups, "=", remaining_dups/dups, "\n", file=outfile, append=TRUE)

    dups <- length(unique(order_clusters$Species_resolved[grepl(" ", order_clusters$Species_resolved) & duplicated(order_clusters$Species_resolved)]))
    correct_species <- order_clusters[grepl(" ", order_clusters$Species_resolved), ]
    non_numt_correct_species <- correct_species[!correct_species$cluster %in% numt_clusters, ]
    remaining_dups <- sum(table(non_numt_correct_species$Species) > 1)
    
    cat ("False negative rate (duplicated resolved species remaining): ", remaining_dups, "of", dups, "=", remaining_dups/dups, "\n", file=outfile, append=TRUE)
  }

  if ("BOLD_bin_resolved" %in% colnames(order_clusters)) {
    unique_bins           <- length(unique(order_clusters$BOLD_bin_resolved[grepl("BOL", order_clusters$BOLD_bin_resolved)]))
    remaining_unique_bins <- length(unique(order_clusters$BOLD_bin_resolved[grepl("BOL", order_clusters$BOLD_bin_resolved) & !(order_clusters$cluster %in% numt_clusters)]))
    cat ("False positive rate (unique resolved BOLD bins removed): ", unique_species-remaining_unique_species, "of", unique_species, "=", (unique_species-remaining_unique_species)/unique_species, "\n", file=outfile, append=TRUE)
    
    dups <- sum(duplicated(order_clusters$BOLD_bin_resolved[grepl("BOL", order_clusters$BOLD_bin_resolved)]))
    remaining_dups <- sum(duplicated(order_clusters$BOLD_bin_resolved[grepl("BOL", order_clusters$BOLD_bin_resolved) & !(order_clusters$cluster %in% numt_clusters)]))
    cat ("False negative rate (duplicated resolved BOLD bin records remaining): ", remaining_dups, "of", dups, "=", remaining_dups/dups, "\n", file=outfile, append=TRUE)
    
    dups <- length(unique(order_clusters$BOLD_bin_resolved[grepl("BOL", order_clusters$BOLD_bin_resolved) & duplicated(order_clusters$BOLD_bin_resolved)]))
    correct_bins <- order_clusters[grepl("BOL", order_clusters$BOLD_bin_resolved), ]
    non_numt_correct_bins <- correct_bins[!correct_bins$cluster %in% numt_clusters, ]
    remaining_dups <- length(unique(non_numt_correct_bins$BOLD_bin_resolved[duplicated(non_numt_correct_bins$BOLD_bin_resolved)]))
    cat ("False negative rate (duplicated resolved BOLD bins remaining): ", remaining_dups, "of", dups, "=", remaining_dups/dups, "\n", file=outfile, append=TRUE)
  }
  
  if ("man_numt" %in% colnames(order_clusters)) {
    man_numt_clusters <- order_clusters$cluster[order_clusters$man_numt]
    n_numts <- length(man_numt_clusters)
    remaining_numts <- n_numts - sum(man_numt_clusters %in% numt_clusters)
    cat ("False negative rate (manually identified numt clusters remaining): ", remaining_numts, "of", n_numts, "=", remaining_numts/n_numts, "\n", file=outfile, append=TRUE)
  }

  if ("Species_updated" %in% colnames(order_clusters)) {
    unique_species <- length(unique(order_clusters$Species_updated[grepl(" ", order_clusters$Species_updated)]))
    remaining_unique_species <- length(unique(order_clusters$Species_updated[grepl(" ", order_clusters$Species_updated) & !(order_clusters$cluster %in% numt_clusters)]))
    cat ("False positive rate (unique updated species removed): ", unique_species-remaining_unique_species, "of", unique_species, "=", (unique_species-remaining_unique_species)/unique_species, "\n", file=outfile, append=TRUE)
    
    dups <- length(unique(order_clusters$Species_updated[grepl(" ", order_clusters$Species_updated) & duplicated(order_clusters$Species_updated)]))
    species_with_space <- order_clusters[grepl(" ", order_clusters$Species_updated), ]
    non_numt_species_with_space <- species_with_space[!species_with_space$cluster %in% numt_clusters, ]
    remaining_dups <- length(unique(non_numt_species_with_space$Species_updated[duplicated(non_numt_species_with_space$Species_updated)]))

    cat ("False negative rate (duplicated updated species remaining): ", remaining_dups, "of", dups, "=", remaining_dups/dups, "\n", file=outfile, append=TRUE)

    if ("man_numt" %in% colnames(order_clusters)) {
      man_numt_clusters <- order_clusters$cluster[order_clusters$man_numt]

      unique_species <- length(unique(order_clusters$Species_updated[grepl(" ", order_clusters$Species_updated) & !(order_clusters$cluster %in% man_numt_clusters)]))
      remaining_unique_species <- length(unique(order_clusters$Species_updated[grepl(" ", order_clusters$Species_updated) & !(order_clusters$cluster %in% man_numt_clusters) & !(order_clusters$cluster %in% numt_clusters)]))
      cat ("False positive rate (unique updated species, cleaned from man numts, removed): ", unique_species-remaining_unique_species, "of", unique_species, "=", (unique_species-remaining_unique_species)/unique_species, "\n", file=outfile, append=TRUE)
      
      dups           <- sum(duplicated(order_clusters$Species_updated[grepl(" ", order_clusters$Species_updated) & !(order_clusters$cluster %in% man_numt_clusters)]))
      remaining_dups <- sum(duplicated(order_clusters$Species_updated[grepl(" ", order_clusters$Species_updated) & !(order_clusters$cluster %in% man_numt_clusters) & !(order_clusters$cluster %in% numt_clusters)]))
      cat ("False negative rate (duplicated updated species records, cleaned from man numts, remaining): ", remaining_dups, "of", dups, "=", remaining_dups/dups, "\n", file=outfile, append=TRUE)
      
      dups           <- length(unique(order_clusters$Species_updated[grepl(" ", order_clusters$Species_updated) & duplicated(order_clusters$Species_updated) & !(order_clusters$cluster %in% man_numt_clusters)]))
      remaining_dups <- length(unique(order_clusters$Species_updated[grepl(" ", order_clusters$Species_updated) & duplicated(order_clusters$Species_updated) & !(order_clusters$cluster %in% man_numt_clusters) & !(order_clusters$cluster %in% numt_clusters)]))
      cat ("False negative rate (duplicated updated species, cleaned from man numts, remaining): ", remaining_dups, "of", dups, "=", remaining_dups/dups, "\n", file=outfile, append=TRUE)
    }
  }
}

list_removed_species <- function(order_clusters, numt_data) {

  numt_clusters <- numt_data$cluster[numt_data$numt]
  man_numt_clusters <- order_clusters$cluster[order_clusters$man_numt]
  
  unique_species <- unique(order_clusters$Species_updated[grepl(" ", order_clusters$Species_updated) & !(order_clusters$cluster %in% man_numt_clusters)])
  remaining_unique_species <- unique(order_clusters$Species_updated[grepl(" ", order_clusters$Species_updated) & !(order_clusters$cluster %in% man_numt_clusters) & !(order_clusters$cluster %in% numt_clusters)])

  removed_species <- unique_species[!(unique_species %in% remaining_unique_species)]
  
  cat ("Removed unique species, cleaned from man numts:\n")
  print(removed_species)

  removed_species
}

complementing_data <- function(order_name, taxfile, consensus_taxfile, countsfile, outfile, threads=1) {
    # Read in nochimera cluster information and create core of the data
    # on the order clusters
    F <- read.table(taxfile,header=TRUE,sep="\t") # taxonomy
    order_clusters <- F[F$Order==order_name & F$representative==1,]
    # Add data on resolved cluster taxonomy
    G <- read.table(consensus_taxfile,header=TRUE,sep="\t")
    order_clusters <- merge(order_clusters, G[, c("cluster","Family","Genus","Species","BOLD_bin")], by="cluster", suffixes = c("","_resolved"), all.x=TRUE)

    # Assume that the original species identification is correct for the representative ASV in the cases where the
    # cluster was filtered out before the name resolution in version 1 of the pipeline. And that the resolved species
    # identification is correct for the remaining ones
    order_clusters$Species_updated <- order_clusters$Species_resolved
    for (i in 1:nrow(order_clusters)) {
      if (is.na(order_clusters$Species_updated[i]))
        order_clusters$Species_updated[i] <- order_clusters$Species[i]
    }

  # Get total number of reads for the clusters
  cat("Reading cluster counts\n")
  D <- fread(countsfile, check.names=FALSE, nThread=threads, 
              sep="\t", header=TRUE, data.table=FALSE, )
  cat(paste0("Subsetting to clusters in", order_name, "\n"))
  rownames(D) <- D$cluster
  D <- D[order_clusters$cluster,2:ncol(D)]

  order_clusters$n_reads <- rowSums(D)
  order_clusters$n_samples <- apply(D[,2:ncol(D)] > 0, 1, sum)
  write.table(order_clusters, outfile, sep="\t", row.names=FALSE, quote=FALSE)
}