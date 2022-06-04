
log <- file(snakemake@log$progress, open="wt")
sink(log)
sink(log, type="message")

if (!"lulu"%in%installed.packages()) {
  library(devtools)
  install_github("tobiasgf/lulu")
}

library(lulu)

tmp <- Sys.getenv("TMPDIR")
tmpdir <- gsub(x=snakemake@params$tmpdir, pattern="^\\$TMPDIR", replacement=tmp)
print(tmpdir)
dist <- gsub(x=snakemake@params$dist, pattern="^\\$TMPDIR", replacement=tmp)

system_string <- paste(c("mkdir -p"), tmpdir)
print(system_string)
system(system_string)

system_string <- paste(c("gunzip -c"), c(snakemake@input$dist), c(">"), dist)
print(system_string)
system(system_string)

otutable <- read.table(snakemake@input$counts, sep="\t", row.names=1, header=TRUE)
matchlist <- read.table(dist, sep="\t")

curated_result <- lulu(otutable=otutable, matchlist=matchlist, minimum_ratio_type=snakemake@params$minimum_ratio_type,
                       minimum_ratio = snakemake@params$minimum_ratio, minimum_match = snakemake@params$minimum_match,
                       minimum_relative_cooccurence = snakemake@params$minimum_relative_cooccurence)

write.table(curated_result$curated_table, snakemake@output$curated_table, sep="\t", quote=FALSE)
write.table(curated_result$otu_map, snakemake@output$otu_map, sep="\t", quote=FALSE)

# Move logfile
system_string <- paste(c("mv lulu.log*"), c(snakemake@log$log))
system(system_string)