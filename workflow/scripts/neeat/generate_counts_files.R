sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)
source(snakemake@params$functions)
library(data.table)

## INPUT
counts_file <- snakemake@input$counts
taxonomy_file <- snakemake@input$taxonomy

## PARAMS
meta_file <- snakemake@params$meta
sample_id_col <- snakemake@params$sample_id_col
sample_type_col <- snakemake@params$sample_type_col
sample_val <- snakemake@params$sample_val

# Read in data
cat("Processing data\n")
#meta <- get_se_meta()

if ((meta_file!=FALSE) && file.exists(meta_file)) {
    meta <- read.delim(meta_file)
    rownames(meta) <- meta[, 1]
    if (sample_id_col != "") {
        rownames(meta) <- meta[, sample_id_col]
    }
    if ((sample_type_col!="") && (length(sample_val)>0)) {
        meta <- meta[meta[,sample_type_col]%in%sample_val,]
    }
    samples <- rownames(meta)
} else {
    meta <- FALSE
    samples <- c()
}

counts <- fread(counts_file, sep="\t")
counts <- as.data.frame(counts)
taxonomy <- read.delim(taxonomy_file)

# Extract counts for clusters in taxonomy
counts <- counts[counts$cluster %in% taxonomy$cluster,]

cat("Generating raw counts table\n")
# extract counts for samples
if (is.null(samples)) {
    index <- 2:ncol(counts)
} else {
    index <- match(samples,colnames(counts))
}
index <- c(1, index[!is.na(index)])
counts <- counts[,index]
# extract counts with rowSum > 0
include_rows <- rowSums(counts[,2:ncol(counts)])>0
counts <- counts[include_rows,]
# extract counts with colSum > 0
include_cols <- c(TRUE,as.logical(colSums(counts[,2:ncol(counts)])>0))
counts <- counts[,include_cols]
#write.table(counts,"../bigdata/cluster_counts.tsv", row.names=FALSE, sep="\t")
write.table(counts, snakemake@output$cluster_counts, row.names=FALSE, sep="\t", quote=FALSE)

quit(save="no", status=0)

sink()