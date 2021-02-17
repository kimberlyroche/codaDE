# Parse "absolute" abundance 16S data from Morton et al. "reference frames" paper

file_dir <- "C:/Users/kim/Documents/codaDE/data/Morton_2019/Github"

otu_table <- read.delim(file.path(file_dir, "oral_trimmed_deblur.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
metadata <- read.delim(file.path(file_dir, "oral_trimmed_metadata.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# "flowcount" is calculated this way in the paper (see Python notebook for Fig. 2)
flowcount <- (metadata$flow.cell.5min.1 + metadata$flow.cell.5min.2) / 2
event_labels <- metadata$brushing_event
sample_labels <- metadata$SampleID

# Compute relative abundances
tax_legend <- otu_table$OTU_ID
props <- otu_table[,2:ncol(otu_table)]
props <- apply(props, 2, function(x) x/sum(x))

# Scale these up to absolute abundances
# `props` is taxa x samples
counts <- props
for(i in 1:ncol(counts)) {
  counts[,i] <- counts[,i]*flowcount[i]
}

# Re-order samples by condition (event before/after)
counts <- counts[,c(which(event_labels == "before"), which(event_labels == "after"))]

# Scale by minimum observed abundance as I did with Barlow et al.
min_observed <- min(counts[which(counts > 0, arr.ind = TRUE)])
counts <- counts / min_observed

saveRDS(list(counts = counts, groups = factor(c(rep("before", sum(event_labels == "before")),
                                                rep("after", sum(event_labels == "after"))), levels = c("before", "after")), tax = NULL),
        file = file.path(file_dir, "absolute_parsed.rds"))
