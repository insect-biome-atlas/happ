sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)
## SOURCE FUNCTIONS
source(snakemake@params$echo_filter)
source(snakemake@params$evo_filter)
source(snakemake@params$abundance_filter)

## WILDCARDS
assignment_rank <- snakemake@wildcards$noise_rank

## PARAMS
min_match <- snakemake@params$min_match
n_closest <- snakemake@params$n_closest
echo_min_overlap <- snakemake@params$echo_min_overlap
echo_read_ratio_type <- snakemake@params$echo_read_ratio_type
echo_max_read_ratio <- snakemake@params$echo_max_read_ratio
echo_require_corr <- snakemake@params$echo_require_corr
evo_local_min_overlap <- snakemake@params$evo_local_min_overlap
dist_type_local <- snakemake@params$dist_type_local
dist_threshold_local <- snakemake@params$dist_threshold_local
dist_threshold_global <- snakemake@params$dist_threshold_global
abundance_cutoff_type <- snakemake@params$abundance_cutoff_type
abundance_cutoff <- snakemake@params$abundance_cutoff
assignment_rank <- snakemake@params$assignment_rank

## INPUT
counts_file <- snakemake@input$counts
dist_file <- snakemake@input$distlist
taxonomy_file <- snakemake@input$taxonomy

counts <- read.delim(counts_file, sep="\t")
rownames(counts) <- counts[, 1]
counts <- counts[, 2:ncol(counts)]
distlist_lines <- length(readLines(dist_file, n=1))
if (distlist_lines == 0) {
     distlist <- data.frame()
} else {
     distlist <- read.delim(dist_file, sep="\t")
}
taxonomy <- read.delim(taxonomy_file, sep="\t")
results <- neeat_filter(
     counts=counts,
     distlist=distlist,
     taxonomy=taxonomy,
     min_match=min_match,
     n_closest=n_closest,
     echo_min_overlap=echo_min_overlap,
     echo_read_ratio_type=echo_read_ratio_type,
     echo_max_read_ratio=echo_max_read_ratio,
     echo_require_corr=echo_require_corr,
     evo_local_min_overlap=evo_local_min_overlap,
     dist_type_local=dist_type_local,
     dist_threshold_local=dist_threshold_local,
     dist_threshold_global=dist_threshold_global,
     abundance_cutoff_type=abundance_cutoff_type,
     abundance_cutoff=abundance_cutoff,
     assignment_rank=assignment_rank
     )

# output retained taxonomy df
rownames(taxonomy) <- taxonomy$ASV
taxonomy_retained <- taxonomy[results$retained_clusters,]
write.table(taxonomy_retained, snakemake@output$retained, sep="\t", row.names=FALSE, quote=FALSE)

# output discarded taxonomy df
taxonomy_discarded <- taxonomy[results$discarded_clusters,]
write.table(taxonomy_discarded, snakemake@output$discarded, sep="\t", row.names=FALSE, quote=FALSE)

# output filtered counts
write.table(results$filtered_counts, snakemake@output$counts, sep="\t", row.names=TRUE, quote=FALSE)

sink()