library(data.table)

source("spikes_controls_fxns.R")

# Generate cleaned filtered counts for Sweden

# Read in filtered counts
cat("Reading in filtered counts for SE\n")
filtered_counts <- fread("../results/filtered_counts_SE.tsv")

# Read in raw counts
cat("Reading in raw cluster counts for SE\n")
raw_counts <- fread("~/dev/figshare-repos/iba/processed_data/SE.v2/cluster_counts.tsv")
idx <- which(colSums(raw_counts[,2:ncol(raw_counts)])!=0)
idx <- c(1,idx+1)
raw_counts <- raw_counts[,..idx]

# Read in taxonomy
taxonomy <- read.delim("~/dev/figshare-repos/iba/processed_data/SE.v2/cluster_taxonomy.tsv")
taxonomy <- taxonomy[taxonomy$representative==1,]

# Read in metadata
meta <- read.delim("~/dev/figshare-repos/iba/raw_data/CO1_sequencing_metadata_SE.tsv")

# Get samples and controls
samples <- meta$sampleID_NGI[meta$lab_sample_type=="sample" & meta$sequencing_status=="sequencing successful"]
control_types <- c("buffer_blank","pcr_neg","extraction_neg","buffer_blank_art_spikes")
controls <- meta$sampleID_NGI[meta$lab_sample_type %in% control_types]
spikein_samples <- meta$sampleID_NGI[meta$dataset %in% c("CO1_lysate_SE","CO1_homogenate_SE")]

cleaned_filtered_counts <- remove_controls(filtered_counts,raw_counts,taxonomy,samples,controls)

cleaned_filtered_counts <- remove_spikes(cleaned_filtered_counts,spikein_samples_taxonomy)

write.table(cleaned_filtered_counts,"../results/cleaned_filtered_counts_SE.tsv")


# Generate cleaned filtered counts for Madagascar

# Read in filtered counts
cat("Reading in filtered counts for MG\n")
filtered_counts <- fread("../results/filtered_counts_MG.tsv")

# Read in raw counts
cat("Reading in raw cluster counts for MG\n")
raw_counts <- fread("~/dev/figshare-repos/iba/processed_data/MG.v2/cluster_counts.corrected.NGI_IDs.tsv")
idx <- which(colSums(raw_counts[,2:ncol(raw_counts)])!=0)
idx <- c(1,idx+1)
raw_counts <- raw_counts[,..idx]

# Read in taxonomy
taxonomy <- read.delim("~/dev/figshare-repos/iba/processed_data/MG.v2/cluster_taxonomy.tsv")
taxonomy <- taxonomy[taxonomy$representative==1,]

# Read in metadata
meta <- read.delim("~/dev/figshare-repos/iba/raw_data/CO1_sequencing_metadata_MG.new.tsv")

# Get samples and controls
samples <- meta$sampleID_NGI[meta$lab_sample_type=="sample" & grepl("successful",meta$sequencing_status)]
control_types <- c("buffer_blank","pcr_neg","extraction_neg")
controls <- meta$sampleID_NGI[meta$lab_sample_type %in% control_types]
spikein_samples <- meta$sampleID_NGI[meta$dataset=="CO1_lysate_MG"]

cleaned_filtered_counts <- remove_controls(filtered_counts,raw_counts,taxonomy,samples,controls)
                     
cleaned_filtered_counts <- remove_spikes(cleaned_filtered_counts,spikein_samples,taxonomy,cutoff=0.80)

write.table(cleaned_filtered_counts,"../results/cleaned_filtered_counts_MG.tsv")

