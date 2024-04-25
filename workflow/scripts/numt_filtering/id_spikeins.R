## ID SPIKEIN CLUSTERS
id_spikein_clusters <- function(taxonomy_file, counts_file, spikein_file, method="species", threads=1) {
    library(data.table)
    # Read in data
    cat ("Loading data\n")
    taxonomy <- read.table(taxonomy_file, header=TRUE, sep="\t")
    colnames(taxonomy) <- tolower(colnames(taxonomy))
    taxonomy <- data.table(taxonomy)
    counts <- fread(counts_file, check.names=FALSE, nThread=threads, 
                sep="\t", header=TRUE, data.table=FALSE)
    spikeins <- read.table(spikein_file, header=TRUE, sep="\t")
    colnames(spikeins) <- tolower(colnames(spikeins))
    cat ("Identifying spikein clusters...\n")
    # ID potential clusters with species name
    if (method == "species") {
        # ID potential clusters with species name
        spikein_tax <- subset(taxonomy, grepl(paste(spikeins$species, collapse='|'), species, ignore.case = TRUE))
    }
    if (method == "BOLD_bin") {
        # ID potential clusters with BOLD bin
        spikein_tax <- subset(taxonomy, grepl(paste(unlist(strsplit(spikeins$bold_bin, "; ")), collapse='|'), bold_bin, ignore.case = TRUE))
    }  

    # Get potential cluster counts
    spikein_counts <- counts[counts$cluster %in% unique(spikein_tax$cluster),]
    
    # Filter potential clusters on presence > 50%
    cluster_presence <- rowMeans(spikein_counts[, -1] > 0)
    spikein_clusters <- spikein_counts$cluster[cluster_presence >= 0.5]   
    spikein_clusters_tax <- taxonomy[cluster%in%spikein_clusters & representative==1]
    print(paste0(length(spikein_clusters), " clusters identified"))
    if ("bold_bin"%in%colnames(spikeins)) {
        uniq_bins <- unique(unlist(strsplit(spikeins$bold, "; ")))
        found_bold_ids <- uniq_bins[uniq_bins%in%spikein_clusters_tax$bold]
        not_found_bold_ids <- setdiff(uniq_bins, found_bold_ids)
        print(paste0(length(found_bold_ids),"/",length(uniq_bins), " BOLD_bin ids found: ", paste0(found_bold_ids, collapse=",")))
        print(paste0(length(not_found_bold_ids), "/", length(uniq_bins), " BOLD_bin ids NOT found: ", paste0(not_found_bold_ids, collapse=",")))
    }
    if ("species"%in%colnames(spikeins)) {
        uniq_species <- unique(spikeins$species)
        found_species <- uniq_species[uniq_species%in%spikein_clusters_tax$species]
        not_found_species <- setdiff(uniq_species, found_species)
        print(paste0(length(found_species), "/", length(uniq_species), " species found: ", paste0(found_species, collapse=",")))
        print(paste0(length(not_found_species), "/", length(uniq_species), " species NOT found: ", paste0(not_found_species, collapse=",")))
    }
    spikein_clusters
 }

sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)

taxonomy_file <- snakemake@input$taxonomy_file
counts_file <- snakemake@input$counts_file
spikein_file <- snakemake@input$spikein_file
method <- snakemake@params$method
threads <- snakemake@threads
outfile <- snakemake@output[[1]]

spikein_clusters <- id_spikein_clusters(taxonomy_file, counts_file, spikein_file, method, threads)
fwrite(list(spikein_clusters), file = opt$outfile)
sink()