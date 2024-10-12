# Abundance filter algorithm
abundance_filter <- function(counts, cutoff, cutoff_type="sum") {

  if (cutoff_type=="sum") {
    reads <- rowSums(counts)
  } else {
    reads <- apply(counts, 1, max)
  }

  retained_clusters <- rownames(counts)[reads>=cutoff]
  discarded_clusters <- rownames(counts)[reads<cutoff]

  list(retained_clusters=retained_clusters, discarded_clusters=discarded_clusters)
}

