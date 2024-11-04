# Evaluate results from cluster (=otu) filtering
# Note that discarded OTUs are assumed to be given as rep ASVs
# and not cluster names 
eval_res <- function(discarded_otus, otu_taxonomy, finbol_tax, asv_tax, se_family, cutoff=0.99, debug=TRUE) {
  
  if (debug) {
    cat ("The taxonomy has", nrow(otu_taxonomy), "clusters in total\n")
    cat (sum(grepl(" ",otu_taxonomy$Species)), "clusters have rep asv annotated to species\n")
    cat ("There are", length(unique(otu_taxonomy$Species[grepl(" ",otu_taxonomy$Species)])), "unique species annotations\n")
    cat (sum(grepl("BOL",otu_taxonomy$BOLD_bin)), "clusters have rep asv annotated to BOLD bin\n")
    cat ("There are", length(unique(otu_taxonomy$BOLD_bin[grepl("BOL",otu_taxonomy$BOLD_bin)])), "unique BOLD bin annotations\n")
    cat ("There are", length(discarded_otus), "clusters removed\n")
  }

  unique_bins <- length(unique(otu_taxonomy$BOLD_bin[grepl("BOL", otu_taxonomy$BOLD_bin)]))
  remaining_unique_bins <- length(unique(otu_taxonomy$BOLD_bin[grepl("BOL", otu_taxonomy$BOLD_bin) & !(otu_taxonomy$ASV %in% discarded_otus)]))
  removed_uniques <- unique_bins - remaining_unique_bins
  false_pos <- removed_uniques / unique_bins
  if (debug) {
    cat ("False positive rate (unique  BOLD bins removed): ", removed_uniques,
         "of", unique_bins, "=", false_pos, "\n")
  }

  n_finbol_clusters <- sum(otu_taxonomy$ASV %in% finbol_tax$rep_asv)
  n_discarded_finbol_clusters <- sum(discarded_otus %in% finbol_tax$rep_asv)
  false_pos_finbol <- n_discarded_finbol_clusters / n_finbol_clusters
  if (debug) {
    cat ("False positive rate (FinBOL clusters removed): ", n_discarded_finbol_clusters,
         "of", n_finbol_clusters, "=", false_pos_finbol, "\n")
  }

  dups <- length(unique(otu_taxonomy$BOLD_bin[grepl("BOL", otu_taxonomy$BOLD_bin) & duplicated(otu_taxonomy$BOLD_bin)]))
  dup_records <- length(otu_taxonomy$BOLD_bin[grepl("BOL", otu_taxonomy$BOLD_bin) & duplicated(otu_taxonomy$BOLD_bin)])
  BOLD_bins <- otu_taxonomy[grepl("BOL", otu_taxonomy$BOLD_bin), ]
  retained_BOLD_bins <- BOLD_bins[!BOLD_bins$ASV %in% discarded_otus, ]
  remaining_dups <- length(unique(retained_BOLD_bins$BOLD_bin[duplicated(retained_BOLD_bins$BOLD_bin)]))
  remaining_dup_records <- sum(duplicated(retained_BOLD_bins$BOLD_bin))
  false_neg = remaining_dups/dups
  if (debug) {
      cat ("False negative rate (duplicated BOLD bins remaining): ", remaining_dups, "of", dups, "=", false_neg, "\n")
  }

  se_known_fams <- se_family$Family[se_family$prop_known > cutoff & se_family$Family %in% unique(otu_taxonomy$Family)]
  known_spp <- sum(se_family$Known_2017[se_family$Family %in% se_known_fams])
  T <- asv_tax[asv_tax$Family %in% se_known_fams & asv_tax$cluster %in% unique(otu_taxonomy$cluster),]
  found_spp <- length(unique(T$Species[grepl(" ",T$Species)]))
  found_spp_clusters <- unique(T$cluster[grepl(" ",T$Species)])
  n_found_spp_clusters <- length(found_spp_clusters)

  clusters_in_known <- unique(otu_taxonomy$cluster[otu_taxonomy$Family %in% se_known_fams])
  n_clusters_in_known <- length(clusters_in_known)
  spurious_clusters <- clusters_in_known[!(clusters_in_known %in% found_spp_clusters)]
  n_spurious <- length(spurious_clusters)
  removed_clusters <- otu_taxonomy$cluster[match(discarded_otus,otu_taxonomy$ASV)]
  removed_spurious <- spurious_clusters[spurious_clusters %in% removed_clusters]
  n_removed_spurious <- length(removed_spurious)
  n_remaining_spurious <- n_spurious - n_removed_spurious
  false_neg_spurious <- n_remaining_spurious / n_spurious
  
  if (debug) {
    cat ("Number of clusters in known: ", n_clusters_in_known, "\n")
    cat ("Number of found spp in known: ", found_spp, "\n")
    cat ("Number of found spp clusters in known: ", n_found_spp_clusters, "\n")
    cat ("Number of spurious clusters: ", n_spurious, "\n")
    cat ("False negative rate (remaining fraction of spurious clusters): ", n_remaining_spurious, "remaining spurious clusters of",
         n_spurious, "spurious clusters =", false_neg_spurious, "fraction remaining spurious clusters\n")
  }

  if (n_finbol_clusters != 0)
    false_pos_comb <- (unique_bins*false_pos + n_finbol_clusters*false_pos_finbol) / (unique_bins + n_finbol_clusters)
  else
    false_pos_comb <- false_pos

  if (n_spurious != 0)
    false_neg_comb <- (dups*false_neg + n_spurious*false_neg_spurious) / (dups + n_spurious)
  else
    false_neg_comb <- false_neg

  list(clusters=nrow(otu_taxonomy),
       removed_clusters=length(discarded_otus),
       unique_bins=unique_bins,
       removed_uniques=removed_uniques,
       duplicate_bins=dups,
       duplicate_bin_records=dup_records,
       remaining_dups=remaining_dups,
       remaining_dup_records=remaining_dup_records,
       false_pos=false_pos,
       false_neg=false_neg,
       finbol_clusters=n_finbol_clusters,
       removed_finbol_clusters=n_discarded_finbol_clusters,
       spp_in_known=known_spp,
       found_spp_in_known=found_spp,
       clusters_in_known=n_clusters_in_known,
       spurious=n_spurious,
       remaining_spurious=n_remaining_spurious,
       false_pos_finbol=false_pos_finbol,
       false_neg_spurious=false_neg_spurious,
       false_pos_comb=false_pos_comb,
       false_neg_comb=false_neg_comb
       )
}


# Summarize results across orders
add_hexapoda_res <- function(res,params) {
 
    # Cycle over runs in params df
    for (i in 1:nrow(params)) {

        x <- res$run==params$run[i]

        # Do the complicated math first
        false_pos <- sum(res[x,"removed_uniques"])/sum(res[x,"unique_bins"])
        false_neg <- sum(res[x,"remaining_dups"])/sum(res[x,"duplicate_bins"])

        false_pos_finbol <- sum(res[x,"removed_finbol_clusters"])/sum(res[x,"finbol_clusters"])
        false_neg_spurious <- sum(res[x,"remaining_spurious"])/sum(res[x,"spurious"])

        # We need not guard against small numbers here in finbol and spurious
        unique_bins <- sum(res[x,"unique_bins"])
        n_finbol_clusters <- sum(res[x,"finbol_clusters"])
        duplicate_bins <- sum(res[x,"duplicate_bins"])
        spurious <- sum(res[x,"spurious"])
        false_pos_comb <- (unique_bins*false_pos + n_finbol_clusters*false_pos_finbol) / (unique_bins + n_finbol_clusters)
        false_neg_comb <- (duplicate_bins*false_neg + spurious*false_neg_spurious) / (duplicate_bins + spurious)

        res <- rbind(res,
                     c(taxon="Hexapoda",
                       as.list(params[i,]),
                       clusters                     = sum(res[x,"clusters"]),
                       removed_clusters             = sum(res[x,"removed_clusters"]),
                       unique_bins                  = sum(res[x,"unique_bins"]),
                       removed_uniques              = sum(res[x,"removed_uniques"]),
                       duplicate_bins               = sum(res[x,"duplicate_bins"]),
                       duplicate_bin_records        = sum(res[x,"duplicate_bin_records"]),
                       remaining_dups               = sum(res[x,"remaining_dups"]),
                       remaining_dup_records        = sum(res[x,"remaining_dup_records"]),
                       false_pos                    = false_pos,
                       false_neg                    = false_neg,
                       finbol_clusters              = sum(res[x,"finbol_clusters"]),
                       removed_finbol_clusters      = sum(res[x,"removed_finbol_clusters"]),
                       spp_in_known                 = sum(res[x,"spp_in_known"]),
                       found_spp_in_known           = sum(res[x,"found_spp_in_known"]),
                       clusters_in_known            = sum(res[x,"clusters_in_known"]),
                       spurious                     = sum(res[x,"spurious"]),
                       remaining_spurious           = sum(res[x,"remaining_spurious"]),
                       false_pos_finbol             = false_pos_finbol,
                       false_neg_spurious           = false_neg_spurious,
                       false_pos_comb               = false_pos_comb,
                       false_neg_comb               = false_neg_comb
                      )
                    )
    }

    res
}


