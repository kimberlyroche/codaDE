source("path_fix.R")

library(codaDE)
library(tidyverse)
library(gridExtra)
library(randomForest)
library(optparse)

option_list = list(
  make_option(c("--dataset"),
              type = "character",
              default = "Barlow",
              help = "data set to use: Muraro, ESCA, Gruen, Athanasiadou, Ferreira, Kimmerling, Hashimshony",
              metavar = "character"),
  make_option(c("--baseline"),
              type = "character",
              default = "self",
              help = "calls to use as a reference: self, oracle",
              metavar = "character"),
  make_option(c("--threshold"),
              type = "numeric",
              default = "5",
              help = "minimum mean abundance to threshold on",
              metavar = "numeric")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

dataset_name <- opt$dataset
use_baseline <- opt$baseline
threshold <- opt$threshold
testing <- FALSE

methods_list <- c("ALDEx2", "DESeq2", "MAST", "scran")

if(!(dataset_name %in% c("Muraro", "ESCA", "Gruen", "Athanasiadou",
                         "Ferreira", "Kimmerling", "Hashimshony"))) {
  stop(paste0("Invalid data set: ", dataset_name, "!\n"))
}

if(!(use_baseline %in% c("self", "oracle"))) {
  stop(paste0("Invalid DA baseline: ", use_baseline, "!\n"))
}

if(threshold < 0) {
  stop(paste0("Invalid threshold: ", threshold, "!\n"))
}

palette <- list(ALDEx2 = "#46A06B",
                DESeq2 = "#FF5733",
                MAST = "#EF82BB",
                NBGLM = "#7E54DE",
                scran = "#E3C012",
                simulated = "#DDDDDD")

# ------------------------------------------------------------------------------
#   Parse and wrangle validation data
# ------------------------------------------------------------------------------

abs_data <- do.call(paste0("parse_", dataset_name), list(absolute = TRUE))
rel_data <- do.call(paste0("parse_", dataset_name), list(absolute = FALSE))

if(testing & nrow(abs_data$counts) > 500) {
  k <- 500
  sample_idx <- sample(1:nrow(abs_data$counts), size = k, replace = FALSE)
  abs_data$counts <- abs_data$counts[sample_idx,]
  rel_data$counts <- rel_data$counts[sample_idx,]
}

# Reorient as (samples x features)
ref_data <- t(abs_data$counts)
data <- t(rel_data$counts)
groups <- abs_data$groups
groups <- factor(groups)

# Subsample if tons of cells/samples
set.seed(100)
pairs <- table(groups)
if(pairs[1] > 100) {
  A_sample <- sample(which(groups == names(pairs)[1]), size = 100)
} else {
  A_sample <- which(groups == names(pairs)[1])
}
if(pairs[2] > 100) {
  B_sample <- sample(which(groups == names(pairs)[2]), size = 100)
} else {
  B_sample <- which(groups == names(pairs)[2])
}
ref_data <- ref_data[c(A_sample, B_sample),]
data <- data[c(A_sample, B_sample),]
groups <- groups[c(A_sample, B_sample)]

# Convert to integers, just for DESeq2
ref_data <- apply(ref_data, c(1,2), as.integer)
data <- apply(data, c(1,2), as.integer)

# Look at sparsity; filter out too-low features
retain_features <- colMeans(ref_data) >= threshold & colMeans(data) >= threshold
ref_data <- ref_data[,retain_features]
data <- data[,retain_features]

cat(paste0("Dataset dimensions: ", nrow(ref_data), " x ", ncol(ref_data), "\n"))
cat(paste0("Percent zeros: ", round(sum(data == 0)/(nrow(data)*ncol(data)), 3)*100, "%\n"))

# ------------------------------------------------------------------------------
#   Call discrepancy by each method; save to file if doesn't already exist
#
#   This can be time-consuming!
# ------------------------------------------------------------------------------

plot_labels <- list(FPR = "specificity (1 - FPR)", TPR = "sensitivity (TPR)")

plot_df <- NULL

for(DE_method in methods_list) {
  # Pull saved calls on this data set x method if these exist
  save_fn <- file.path("output",
                       "predictive_fits",
                       "all",
                       paste0("calls_",
                              use_baseline,
                              "_",
                              DE_method,
                              "_",
                              dataset_name,
                              "_threshold",
                              threshold,
                              ".rds"))
  if(!file.exists(save_fn)) {
    cat(paste0("Making calls for ", DE_method, " on ", dataset_name, "...\n"))
    # Get baseline differential abundance calls
    if(use_baseline == "oracle") {
      oracle_calls <- call_DA_NB(ref_data, groups)$pval
    } else {
      oracle_calls <- NULL
    }
    all_calls <- DA_wrapper(ref_data, data, groups, DE_method, oracle_calls)
    rates <- calc_DA_discrepancy(all_calls$calls, all_calls$oracle_calls)
    save_obj <- list(all_calls = all_calls, rates = rates)
    saveRDS(save_obj, save_fn)
  } else {
    score_obj <- readRDS(save_fn)
    all_calls <- score_obj$all_calls
    rates <- score_obj$rates
    plot_df <- rbind(plot_df,
                     data.frame(method = DE_method,
                                TPR = rates$TPR,
                                FPR = rates$FPR))
  }
}

pl <- ggplot(plot_df, aes(x = 1 - FPR, y = TPR, fill = method)) +
  geom_point(size = 3, shape = 21) +
  scale_fill_manual(values = palette) +
  theme_bw() +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  labs(x = paste0("observed ", plot_labels$FPR),
       y = paste0("predicted ", plot_labels$TPR),
       fill = "Data type") +
  theme(legend.position = "none")
ggsave(file.path("output",
                 "images",
                 paste0("validations_",
                        use_baseline,
                        "_",
                        dataset_name,
                        "_threshold",
                        threshold,
                        ".png")),
       plot = pl,
       dpi = 100,
       units = "in",
       height = 4,
       width = 4)
