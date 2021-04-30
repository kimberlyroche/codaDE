library(tidyverse)
library(codaDE)

wrangle_Barlow_data <- function(data) {
  # Clean up data
  md_labels <- c("Diet", "Site", "Day", "mouse", "Cage")
  metadata <- data[,md_labels]
  counts <- data[,!(colnames(data) %in% md_labels)]
  counts <- counts[,2:ncol(counts)]
  tax <- colnames(counts)
  counts <- as.matrix(counts)
  colnames(counts) <- NULL
  rownames(counts) <- NULL
  return(list(counts = counts, metadata = metadata, tax = tax))
}

# ------------------------------------------------------------------------------
#   Parse Barlow et al. data
# ------------------------------------------------------------------------------

file_dir <- file.path("data", "Barlow_2020")
absolute_data <- read.table(file.path(file_dir, "Absolute_Abundance_Table.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
relative_data <- read.table(file.path(file_dir, "Relative_Abundance_Table.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)

absolute_obj <- wrangle_Barlow_data(absolute_data)
absolute_counts <- absolute_obj$counts

relative_obj <- wrangle_Barlow_data(relative_data)
relative_counts <- relative_obj$counts

if(any(absolute_obj$tax != relative_obj$tax)) {
  stop("Taxonomy doesn't match between absolute and relative abundances!")
}
tax <- absolute_obj$tax

if(any(absolute_obj$metadata != relative_obj$metadata)) {
  stop("Metadata doesn't match between absolute and relative abundances!")
}
metadata <- absolute_obj$metadata

# Note relative_counts are scaled to 100; scale to CPM
relative_counts <- relative_counts * (1e6 / sum(relative_counts[1,]))

# Pull stool samples from days 4, 7, 10
keto_idx <- which(metadata$Diet == "Keto" & metadata$Day %in% c(4,7,10) & metadata$Site == "Stool")
ctrl_idx <- which(metadata$Diet == "Control" & metadata$Day %in% c(4,7,10) & metadata$Site == "Stool")

days <- metadata[c(ctrl_idx,keto_idx),]$Day

absolute_counts <- rbind(absolute_counts[ctrl_idx,], absolute_counts[keto_idx,])
relative_counts <- rbind(relative_counts[ctrl_idx,], relative_counts[keto_idx,])
labels <- c(rep("control", length(ctrl_idx)), rep("keto", length(keto_idx)))
groups <- labels
groups[groups == "control"] <- 0
groups[groups == "keto"] <- 1
groups <- as.numeric(groups)

# Remove unobserved taxa
absolute_missing <- colSums(absolute_counts) == 0
relative_missing <- colSums(relative_counts) == 0
if(any(absolute_missing != relative_missing)) {
  stop("Missing taxa in absolute vs. relative abundances don't agree!")
}
absolute_counts <- absolute_counts[,!absolute_missing]
relative_counts <- relative_counts[,!relative_missing]
tax <- tax[!absolute_missing]

# Scale each set of counts such that minimum observed count is 1
min_observed <- min(absolute_counts[absolute_counts != 0])
absolute_counts <- absolute_counts/min_observed
min_observed <- min(relative_counts[relative_counts != 0])
relative_counts <- relative_counts/min_observed

# Round all; some methods don't tolerate floats
absolute_counts <- ceiling(absolute_counts)
relative_counts <- ceiling(relative_counts)

# Wrangle data for plotting
absolute_df <- cbind(1:nrow(absolute_counts), absolute_counts)
colnames(absolute_df) <- c("sample", 1:length(tax))
absolute_df <- pivot_longer(as.data.frame(absolute_df), !sample, names_to = "taxon", values_to = "counts")
ggplot(absolute_df, aes(x = log(counts + 1))) +
  geom_histogram()

relative_df <- cbind(1:nrow(relative_counts), relative_counts)
colnames(relative_df) <- c("sample", 1:length(tax))
relative_df <- pivot_longer(as.data.frame(relative_df), !sample, names_to = "taxon", values_to = "relative_counts")
ggplot(relative_df, aes(x = log(relative_counts + 1))) +
  geom_histogram()

# Histograms between conditions
plot_df <- data.frame(count = rowSums(absolute_counts[groups == 0,]),
                      sample = paste0("control", 1:sum(groups == 0)),
                      condition = "control diet")
plot_df <- rbind(plot_df,
                 data.frame(count = rowSums(absolute_counts[groups == 1,]),
                            sample = paste0("keto", 1:sum(groups == 1)),
                            condition = "ketogenic diet"))
plot_df
ggplot(plot_df, aes(x = sample, y = count, fill = factor(condition))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#43C59E", "#7B7CA7")) +
  labs(x = "sample index", y = "reconstructed total abundance", fill = "condition") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave(file.path("output",
                 "images",
                 paste0("Barlow_totals.png")),
       units = "in",
       height = 3,
       width = 5)

# rows are samples (103), columns are taxa (141)

# ------------------------------------------------------------------------------
#   Apply NB GLM to absolute counts to get baseline DE
# ------------------------------------------------------------------------------

oracle_calls <- sapply(1:ncol(absolute_counts), function(idx) { call_DA_NB(absolute_counts, groups, idx) } )
oracle_calls <- oracle_calls < 0.05

call_DA_validation <- function(absolute_counts) {
  # NB GLM
  NB_calls <- sapply(1:ncol(absolute_counts), function(idx) { call_DA_NB(relative_counts, groups, idx) } )
  # DESeq2
  DESeq2_calls <- call_DA_Seurat(relative_counts, groups, method = "DESeq2")
  # MAST
  MAST_calls <- call_DA_MAST(relative_counts, groups)
  # scran
  scran_calls <- call_DA_scran(relative_counts, groups)
  # ALDEx2
  ALDEx2_calls <- call_DA_ALDEx2(relative_counts, groups)
  return(list(baseline = NB_calls,
              DESeq2 = DESeq2_calls,
              MAST = MAST_calls,
              scran = scran_calls,
              ALDEx2 = ALDEx2_calls))
}

calls_list <- call_DA_validation(absolute_counts)

calculate_rates_validation <- function(calls_list) {
  results <- data.frame(tpr = c(),
                        fpr = c(),
                        method = c())
  for(method_idx in 1:length(calls_list)) {
    DE_calls <- calls_list[[method_idx]]
    DE_calls <- DE_calls < 0.05
    
    TP <- sum(oracle_calls & DE_calls)
    FP <- sum(!oracle_calls & DE_calls)
    
    TN <- sum(!oracle_calls & !DE_calls)
    FN <- sum(oracle_calls & !DE_calls)
    
    TPR <- TP / (TP + FN)
    FPR <- FP / (FP + TN)
    
    results <- rbind(results,
                     data.frame(tpr = TPR,
                                fpr = FPR,
                                method = names(calls_list)[method_idx]))
  }
  return(results)
}

results <- calculate_rates_validation(calls_list)
saveRDS(results, file = file.path("output", "Barlow_validation_results.rds"))

# This can be plotted using plot_compiled.R
