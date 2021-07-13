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
              help = "data set to use: Barlow, Morton, Athanasiadou_yeast, Athanasiadou_ciona, TCGA_ESCA",
              metavar = "character"),
  make_option(c("--baseline"),
              type = "character",
              default = "self",
              help = "calls to use as a reference: self, threshold",
              metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

dataset_name <- opt$dataset
use_baseline <- opt$baseline

if(!(dataset_name %in% c("Barlow", "Morton", "Athanasiadou_yeast",
                         "Athanasiadou_ciona", "TCGA_ESCA"))) {
  stop(paste0("Invalid data set: ", dataset_name, "!\n"))
}

if(!(use_baseline %in% c("self", "threshold"))) {
  stop(paste0("Invalid DA baseline: ", use_baseline, "!\n"))
}

palette <- list(ALDEx2 = "#46A06B",
                DESeq2 = "#FF5733",
                MAST = "#EF82BB",
                NBGLM = "#7E54DE",
                scran = "#E3C012",
                simulated = "#DDDDDD")

model <- "RF"

# ------------------------------------------------------------------------------
#   Parse and wrangle validation data
# ------------------------------------------------------------------------------

if(dataset_name == "Barlow") {
  abs_data <- parse_Barlow(absolute = TRUE)
  rel_data <- parse_Barlow(absolute = FALSE)
}
if(dataset_name == "Morton") {
  abs_data <- parse_Morton(absolute = TRUE)
  rel_data <- parse_Morton(absolute = FALSE)
}
if(dataset_name == "Athanasiadou_ciona") {
  abs_data <- parse_Athanasiadou(absolute = TRUE, which_data = "ciona")
  rel_data <- parse_Athanasiadou(absolute = FALSE, which_data = "ciona")
  # Downsample for testing
  # k <- 1000
  # sample_idx <- sample(1:nrow(abs_data$counts), size = k, replace = FALSE)
  # abs_data$counts <- abs_data$counts[sample_idx,]
  # rel_data$counts <- rel_data$counts[sample_idx,]
}
if(dataset_name == "Athanasiadou_yeast") {
  abs_data <- parse_Athanasiadou(absolute = TRUE, which_data = "yeast")
  rel_data <- parse_Athanasiadou(absolute = FALSE, which_data = "yeast")
  # Downsample for testing
  # k <- 1000
  # sample_idx <- sample(1:nrow(abs_data$counts), size = k, replace = FALSE)
  # abs_data$counts <- abs_data$counts[sample_idx,]
  # rel_data$counts <- rel_data$counts[sample_idx,]
}
if(dataset_name == "TCGA_ESCA") {
  abs_data <- parse_ESCA(absolute = TRUE)
  rel_data <- parse_ESCA(absolute = FALSE)
  # Absolute data is hugely lower in abundance, scale it up for testing
  abs_data$counts <- abs_data$counts * 1e05
  # Downsample for testing
  k <- 1000
  sample_idx <- sample(1:nrow(abs_data$counts), size = k, replace = FALSE)
  abs_data$counts <- abs_data$counts[sample_idx,]
  rel_data$counts <- rel_data$counts[sample_idx,]
}

# Reorient as (samples x features)
ref_data <- t(abs_data$counts)
data <- t(rel_data$counts)
groups <- abs_data$groups

if(dataset_name == "Athanasiadou_ciona") {
  ref_data <- ref_data[groups %in% c("lacz", "dnfgfr"),]
  data <- data[groups %in% c("lacz", "dnfgfr"),]
  groups <- factor(groups[groups %in% c("lacz", "dnfgfr")])
}
if(dataset_name == "Athanasiadou_yeast") {
  ref_data <- ref_data[groups %in% c("C12", "C30"),]
  data <- data[groups %in% c("C12", "C30"),]
  groups <- factor(groups[groups %in% c("C12", "C30")])
}

# Convert to integers, just for DESeq2
ref_data <- apply(ref_data, c(1,2), as.integer)
data <- apply(data, c(1,2), as.integer)

# Look at sparsity; filter out too-low features
retain_features <- colSums(ref_data) > 10 & colSums(data) > 10
ref_data <- ref_data[,retain_features]
data <- data[,retain_features]

cat(paste0("Dataset dimensions: ", nrow(ref_data), " x ", ncol(ref_data), "\n"))
cat(paste0("Percent zeros: ", round(sum(data == 0)/(nrow(data)*ncol(data)), 3)*100, "%\n"))

# ------------------------------------------------------------------------------
#   Wrangle data for prediction-making
# ------------------------------------------------------------------------------

# Get info we need from data to make a prediction
if(dataset_name == "Barlow") {
  counts_A <- data[groups == "control",]
  counts_B <- data[groups == "keto",]
}
if(dataset_name == "Morton") {
  counts_A <- data[groups == "before",]
  counts_B <- data[groups == "after",]
}
if(dataset_name == "Athanasiadou_ciona") {
  counts_A <- data[groups == "lacz",]
  counts_B <- data[groups == "dnfgfr",]
}
if(dataset_name == "Athanasiadou_yeast") {
  counts_A <- data[groups == "C12",]
  counts_B <- data[groups == "C30",]
}
if(dataset_name == "TCGA_ESCA") {
  counts_A <- data[groups == "SCC",]
  counts_B <- data[groups == "Other",]
}

# Spike-in a minimal observation for any all-zero features
# min_observed <- min(c(min(c(counts_A[counts_A != 0])),
#                       min(c(counts_B[counts_B != 0]))))
# spike_idx <- which(rowSums(counts_A) == 0)
# for(idx in spike_idx) {
#   counts_A[idx,sample(1:ncol(counts_A), size = 1)] <- min_observed
# }
# spike_idx <- which(rowSums(counts_B) == 0)
# for(idx in spike_idx) {
#   counts_B[idx,sample(1:ncol(counts_B), size = 1)] <- min_observed
# }

# This takes 2-3 min. to run on 15K features
features <- characterize_dataset(counts_A, counts_B)

# Add a few more
features$P <- ncol(counts_A)

# ------------------------------------------------------------------------------
#   Make predictions on simulations and real data and visualize these together
# ------------------------------------------------------------------------------

plot_labels <- list(fpr = "specificity (1 - FPR)", tpr = "sensitivity (TPR)")

for(use_result_type in c("tpr", "fpr")) {

  plot_df <- NULL
  
  for(DE_method in c("ALDEx2", "DESeq2", "MAST", "NBGLM", "scran")) {
    
    if(max(table(groups)) < 5 && DE_method == "scran") {
      next
    }
    
    # --------------------------------------------------------------------------
    #   Call discrepancy by chosen method
    # --------------------------------------------------------------------------
    
    # Get baseline differential abundance calls
    if(use_baseline == "threshold") {
      # "Differential" features will be those with mean fold change in abundance of 
      # >= 1.5 or <= 0.5
      oracle_calls <- calc_threshold_DA(ref_data,
                                        nA = sum(groups == groups[1]))
    } else {
      oracle_calls <- NULL
    }
    
    rates <- calc_DE_discrepancy(ref_data,
                                 data,
                                 groups,
                                 method = DE_method,
                                 oracle_calls = oracle_calls)
    
    
    rates
    
    # De-noise the absolute data
    # new_totals <- c(round(rnorm(sum(groups == "before"), mean(rowSums(ref_data)[groups == "before"]), 10000)),
    #                 round(rnorm(sum(groups == "after"), mean(rowSums(ref_data)[groups == "after"]), 10000)))
    # new_ref <- sapply(1:nrow(ref_data), function(x) {
    #   y <- ref_data[x,]
    #   (y/sum(y))*new_totals[x]
    # })
    # new_ref <- t(new_ref)
    # new_ref <- apply(new_ref, c(1,2), as.integer)
    new_ref <- ref_data
    
    # Increase noise in the relative data
    # new_totals <- sapply(sample(rowSums(data)), function(x) min(abs(x-10000), 10000))
    # new_totals <- round(runif(nrow(data), min = 1000, max = 2500))
    # new_alt <- sapply(1:nrow(data), function(x) {
    #   y <- data[x,]
    #   (y/sum(y))*new_totals[x]
    # })
    # new_alt <- t(new_alt)
    # new_alt <- apply(new_alt, c(1,2), as.integer)
    new_alt <- data
    
    # What do observed / reconstructed totals look like?
    plot_df <- data.frame(x = 1:nrow(new_ref),
                          y = rowSums(new_ref),
                          group = groups,
                          type = "absolute")
    plot_df <- rbind(plot_df,
                     data.frame(x = 1:nrow(new_alt),
                                y = rowSums(new_alt),
                                group = groups,
                                type = "relative"))
    ggplot(plot_df, aes(x = x, y = y, fill = group)) +
      geom_bar(stat = "identity") +
      facet_wrap(. ~ type, scales = "free_y") +
      labs(x = "sample index",
           y = "abundance") +
      theme_bw()
    
    # oracle_calls <- calc_threshold_DA(ref_data,
    #                                   nA = sum(groups == groups[1]))
    # oracle_calls <- DA_by_DESeq2(new_ref, data, groups, oracle_calls = NULL)$oracle_calls$pval
    oracle_calls <- NULL
    DE_calls <- DA_by_MAST(new_ref, data, groups, oracle_calls = oracle_calls)
    if(!is.null(oracle_calls)) {
      oracle <- ifelse(oracle_calls < 0.05, 0, 1)
    } else {
      oracle <- ifelse(DE_calls$oracle_calls$pval < 0.05, 0, 1)
    }
    calls <- ifelse(DE_calls$calls$pval < 0.05, 0, 1)

    TP <- sum(oracle == 0 & calls == 0)
    FP <- sum(oracle == 1 & calls == 0)
    TN <- sum(oracle == 1 & calls == 1)
    FN <- sum(oracle == 0 & calls == 1)
    
    TP / (TP + FN)
    
    idx <- which(oracle == 1 & calls == 0)
    
    use_idx <- sample(idx, size = 1)
    plot_df <- data.frame(x = 1:nrow(new_ref),
                          y = new_ref[,use_idx],
                          group = groups,
                          type = "absolute")
    plot_df <- rbind(plot_df,
                     data.frame(x = 1:nrow(data),
                                y = data[,use_idx],
                                group = groups,
                                type = "relative"))
    ggplot(plot_df, aes(x = x, y = y, fill = group)) +
      geom_point(size = 3, shape = 21) +
      facet_wrap(. ~ type, scales = "free_y") +
      labs(x = "sample index",
           y = "abundance") +
      theme_bw()
    
    
    

    # --------------------------------------------------------------------------
    #   Iterate TPR, FPR predictions
    # --------------------------------------------------------------------------
    
    model_fn <- file.path("output",
                          "images",
                          paste0(model, "_results"),
                          "all",
                          paste0(model,
                                 "_all",
                                 "_",
                                 use_result_type,
                                 "_",
                                 use_baseline,
                                 ".rds"))

    if(!file.exists(model_fn)) {
      stop(paste0("Predictive model fit not found: ", model_fn, "!\n"))
    }
    fit_obj <- readRDS(model_fn)
    
    # --------------------------------------------------------------------------
    #   Finish feature wrangling
    # --------------------------------------------------------------------------
    
    features$CORRP <- 1 # Correlated features? (probably)
    features$CORRP <- factor(features$CORRP,
                             levels = levels(fit_obj$train_features$CORRP))
    features$PARTIAL <- 0 # No partial
    features$PARTIAL <- factor(features$PARTIAL,
                               levels = levels(fit_obj$train_features$PARTIAL))
    features$METHOD <- DE_method
    features$METHOD <- factor(features$METHOD,
                              levels = c("ALDEx2",
                                         "DESeq2",
                                         "MAST",
                                         "NBGLM",
                                         "scran"))
    
    # Convert to data.frame
    features_df <- as.data.frame(features)
    
    # Remove some features we're no longer using (too correlated with others)
    features_df <- features_df %>%
      select(-c(FW_RA_PFC1_D, FW_CLR_MED_D, FW_CLR_SD_D, FW_CLR_PNEG_D))
    
    # --------------------------------------------------------------------------
    #   Make predictions on simulated and real
    # --------------------------------------------------------------------------
    
    # if(is.null(plot_df)) {
    #   pred_sim <- predict(fit_obj$result, newdata = fit_obj$test_features)
    #   
    #   plot_df <- data.frame(true = fit_obj$test_response,
    #                         predicted = pred_sim,
    #                         type = "simulated")
    # }

    pred_real <- predict(fit_obj$result, newdata = features_df)
    
    plot_df <- rbind(plot_df,
                     data.frame(true = ifelse(use_result_type == "tpr", rates$tpr, 1 - rates$fpr),
                                predicted = pred_real,
                                type = DE_method))
    
  }
  
  pl <- ggplot() +
    geom_segment(data = data.frame(x = 0, xend = 1, y = 0, yend = 1),
                 mapping = aes(x = x, xend = xend, y = y, yend = yend)) +
    geom_point(data = plot_df,
               mapping = aes(x = true, y = predicted, fill = type),
               shape = 21,
               size = 3) +
    scale_fill_manual(values = palette) + 
    xlim(c(0,1)) +
    ylim(c(0,1)) +
    labs(x = paste0("observed ", plot_labels[[use_result_type]]),
         y = paste0("predicted ", plot_labels[[use_result_type]]),
         fill = "Data type")
  show(pl)
  ggsave(file.path("output",
                   "images",
                   paste0("validations_",
                          use_result_type,
                          "_",
                          use_baseline,
                          "-",
                          dataset_name,
                          ".png")),
         plot = pl,
         dpi = 100,
         units = "in",
         height = 4,
         width = 5.5)
}
