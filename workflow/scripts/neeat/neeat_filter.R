sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)
# NEEAT filter algorithm
# Note that taxonomy needs to be a data frame
# and not a data table for Step 5 code to work
neeat_filter <- function(counts,
                         distlist,
                         taxonomy,
                         steps="echo|evo_local|evo_global|abundance|taxonomy",
                         min_match=84,
                         n_closest=10,
                         echo_min_overlap=0.90,
                         echo_read_ratio_type="mean",
                         echo_max_read_ratio=0.1,
                         echo_require_corr=FALSE,
                         evo_local_min_overlap=0.95,
                         dist_type_local="dadn",
                         dist_threshold_local=1.8,
                         dist_threshold_global=4.0,
                         abundance_cutoff_type="sum",
                         abundance_cutoff=4,
                         assignment_taxonomy=taxonomy,
                         assignment_rank="Order",
                         debug=FALSE,
                         dump_prefix=""
                        ) {

    discarded_clusters <- character()
    all_clusters <- row.names(counts)
    retained_clusters <- all_clusters

    # Step 1. The echo filter
    cat("Running echo filter\n")
    if (grepl("echo",steps) && nrow(counts)>=2 && nrow(distlist)!=0) {
        res <- echo_filter(counts,
                           distlist,
                           min_match=min_match,
                           n_closest=n_closest,
                           min_overlap=echo_min_overlap,
                           read_ratio_type=echo_read_ratio_type,
                           max_read_ratio=echo_max_read_ratio,
                           require_corr=echo_require_corr)

        discarded_clusters <- res$discarded_clusters
        retained_clusters <- res$retained_clusters
        counts <- counts[rownames(counts) %in% retained_clusters,]
    }

    if (debug && length(discarded_clusters)>0) {
        file = paste(dump_prefix,"step1_discarded_otus.tsv",sep="_")
        cat("Writing to file",file,"\n")
        cat(discarded_clusters,sep="\n",file=file,append=TRUE)
    }
    
    cat("Running evo_local filter\n")
    # Step 2. The evo_local filter
    if (grepl("evo_local",steps) && nrow(counts)>=2 && nrow(distlist)!=0) {
        res <- evo_filter(counts,
                          distlist,
                          min_match=min_match,
                          n_closest=n_closest,
                          min_overlap=evo_local_min_overlap,
                          dist_type=dist_type_local,
                          dist_threshold=dist_threshold_local,
                          require_overlap=TRUE)
    
        discarded_clusters <- c(discarded_clusters,res$discarded_clusters)
        retained_clusters <- res$retained_clusters
        counts <- counts[rownames(counts) %in% retained_clusters,]
    }

    if (debug && length(discarded_clusters)>0) {
        file = paste(dump_prefix,"step2_discarded_otus.tsv",sep="_")
        cat("Writing to file",file,"\n")
        cat(discarded_clusters,sep="\n",file=file,append=TRUE)
    }

    # Step 3. The evo_global filter
    if (grepl("evo_global",steps) && nrow(counts)>=2 && nrow(distlist)!=0) {
        res <- evo_filter(counts,
                          distlist,
                          min_match=min_match,
                          n_closest=n_closest,
                          dist_type="wdadn",
                          dist_threshold=dist_threshold_global,
                          require_overlap=FALSE)
    
        discarded_clusters <- c(discarded_clusters,res$discarded_clusters)
        retained_clusters <- res$retained_clusters
        counts <- counts[rownames(counts) %in% retained_clusters,]
    }

    if (debug && length(discarded_clusters)>0) {
        file = paste(dump_prefix,"step3_discarded_otus.tsv",sep="_")
        cat("Writing to file",file,"\n")
        cat(discarded_clusters,sep="\n",file=file,append=TRUE)
    }

    # Step 4: The abundance filter
    if (grepl("abundance",steps)) {
        res <- abundance_filter(counts,
                                cutoff=abundance_cutoff,
                                cutoff_type=abundance_cutoff_type)
                        
        discarded_clusters <- c(discarded_clusters,res$discarded_clusters)
        retained_clusters <- res$retained_clusters
        counts <- counts[rownames(counts) %in% retained_clusters,]
    }
    
    if (debug && length(discarded_clusters)>0) {
        file = paste(dump_prefix,"step4_discarded_otus.tsv",sep="_")
        cat("Writing to file",file,"\n")
        cat(discarded_clusters,sep="\n",file=file,append=TRUE)
    }

    # Step 5: The taxonomy filter
    if (grepl("taxonomy",steps)) {

        unclassified_clusters <- assignment_taxonomy$ASV[grepl("_X",assignment_taxonomy[,assignment_rank]) | grepl("unclassified",assignment_taxonomy[,assignment_rank])]

        retained_clusters <- retained_clusters[!(retained_clusters %in% unclassified_clusters)]
        discarded_clusters <- c(discarded_clusters,retained_clusters[retained_clusters %in% unclassified_clusters])
        counts <- counts[rownames(counts) %in% retained_clusters,]
    }

    if (debug && length(discarded_clusters)>0) {
        file = paste(dump_prefix,"step5_discarded_otus.tsv",sep="_")
        cat("Writing to file",file,"\n")
        cat(discarded_clusters,sep="\n",file=file,append=TRUE)
    }

    # Return retained and discarded clusters, and filtered counts
    list(retained_clusters=retained_clusters, discarded_clusters=discarded_clusters, filtered_counts=counts)
}

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
steps="echo|evo_local|evo_global|abundance|taxonomy",
                         min_match=84,
                         n_closest=10,
                         echo_min_overlap=0.90,
                         echo_read_ratio_type="mean",
                         echo_max_read_ratio=0.1,
                         echo_require_corr=FALSE,
                         evo_local_min_overlap=0.95,
                         dist_type_local="dadn",
                         dist_threshold_local=1.8,
                         dist_threshold_global=4.0,
                         abundance_cutoff_type="sum",
                         abundance_cutoff=4,
                         assignment_taxonomy=taxonomy,
                         assignment_rank="Order",
                         debug=FALSE,
                         dump_prefix=""

## INPUT
counts_file <- snakemake@input$counts
dist_file <- snakemake@input$distlist
taxonomy_file <- snakemake@input$taxonomy

counts <- read.delim(counts_file, sep="\t")
distlist <- read.delim(dist_file, sep="\t")
taxonomy <- read.delim(taxonomy_file, sep="\t")
results <- neeat_filter(counts=counts,
                        distlist=distlist,
                        taxonomy=taxonomy,)
sink()