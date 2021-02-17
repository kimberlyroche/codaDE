# Parse "absolute" abundance 16S data from Barlow et al. paper

file_dir <- "C:/Users/kim/Documents/codaDE/data/Barlow_2020"

data <- read.table(file.path(file_dir, "Absolute_Abundance_Table.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Clean up data
md_labels <- c("Diet", "Site", "Day", "mouse", "Cage")
metadata <- data[,md_labels]
counts <- data[,!(colnames(data) %in% md_labels)]
counts <- counts[,2:ncol(counts)]
tax <- colnames(counts)
counts <- as.matrix(counts)
colnames(counts) <- NULL
rownames(counts) <- NULL

# `counts` is initially 103 samples x 142 taxa

# Pull stool samples from days 4, 7, 10
keto_idx <- which(metadata$Diet == "Keto" & metadata$Day %in% c(4,7,10) & metadata$Site == "Stool")
ctrl_idx <- which(metadata$Diet == "Control" & metadata$Day %in% c(4,7,10) & metadata$Site == "Stool")

counts <- rbind(counts[ctrl_idx,], counts[keto_idx,])
labels <- c(rep("control", length(ctrl_idx)), rep("keto", length(keto_idx)))

# Eliminate all-zero taxa
absent_tax_idx <- which(colSums(counts) == 0)
counts <- counts[,-absent_tax_idx]
tax <- tax[-absent_tax_idx]

# The scaling does weird things to this data set. The minimum *observed* (non-zero) abundance is huge
# giving a gap between observed and unobserved features that is enormous.
# I'm scaling down all the counts by this minimum observed abundance for now. (2/16/2021)
min_observed <- min(counts[which(counts > 0, arr.ind = TRUE)])
counts <- counts / min_observed

saveRDS(list(counts = t(counts), groups = factor(labels, levels = c("control", "keto")), tax = tax), file = file.path(file_dir, "absolute_parsed.rds"))
