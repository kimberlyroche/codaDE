source("path_fix.R")

library(codaDE)
library(tidyverse)
library(gridExtra)
library(randomForest)

palette <- list(ALDEx2 = "#46A06B",
                DESeq2 = "#FF5733",
                MAST = "#EF82BB",
                NBGLM = "#7E54DE",
                scran = "#E3C012",
                simulated = "#DDDDDD")

model <- "RF"
use_baseline <- "self"

# ------------------------------------------------------------------------------
#   Parse and wrangle validation data
# ------------------------------------------------------------------------------

# dataset_name <- "Barlow"
# dataset_name <- "Athanasiadou_ciona"
dataset_name <- "Athanasiadou_yeast"

if(dataset_name == "Barlow") {
  abs_data <- parse_Barlow(absolute = TRUE)
  rel_data <- parse_Barlow(absolute = FALSE)
}
if(dataset_name == "Athanasiadou_ciona") {
  abs_data <- parse_Athanasiadou(absolute = TRUE, which_data = "ciona")
  rel_data <- parse_Athanasiadou(absolute = FALSE, which_data = "ciona")
  # Downsample for testing
  k <- 1000
  sample_idx <- sample(1:nrow(abs_data$counts), size = k, replace = FALSE)
  abs_data$counts <- abs_data$counts[sample_idx,]
  rel_data$counts <- rel_data$counts[sample_idx,]
}
if(dataset_name == "Athanasiadou_yeast") {
  abs_data <- parse_Athanasiadou(absolute = TRUE, which_data = "yeast")
  rel_data <- parse_Athanasiadou(absolute = FALSE, which_data = "yeast")
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

# ------------------------------------------------------------------------------
#   Wrangle data for prediction-making
# ------------------------------------------------------------------------------

# Get info we need from data to make a prediction
if(dataset_name == "Barlow") {
  counts_A <- rel_data$counts[,rel_data$groups == "control"]
  counts_B <- rel_data$counts[,rel_data$groups == "keto"]
}
if(dataset_name == "Athanasiadou_ciona") {
  counts_A <- rel_data$counts[,rel_data$groups == "lacz"]
  counts_B <- rel_data$counts[,rel_data$groups == "dnfgfr"]
}
if(dataset_name == "Athanasiadou_yeast") {
  counts_A <- rel_data$counts[,rel_data$groups == "C12"]
  counts_B <- rel_data$counts[,rel_data$groups == "C30"]
}

# Spike-in a minimal observation for any all-zero features
min_observed <- min(c(min(c(counts_A[counts_A != 0])),
                      min(c(counts_B[counts_B != 0]))))
spike_idx <- which(rowSums(counts_A) == 0)
for(idx in spike_idx) {
  counts_A[idx,sample(1:ncol(counts_A), size = 1)] <- min_observed
}
spike_idx <- which(rowSums(counts_B) == 0)
for(idx in spike_idx) {
  counts_B[idx,sample(1:ncol(counts_B), size = 1)] <- min_observed
}

# This takes 2-3 min. to run on 15K features
features <- characterize_dataset(t(counts_A), t(counts_B))

# Add a few more
features$P <- nrow(counts_A)

# ------------------------------------------------------------------------------
#   Make predictions on simulations and real data and visualize these together
# ------------------------------------------------------------------------------

plot_labels <- list(fpr = "specificity (1 - FPR)", tpr = "sensitivity (TPR)")

for(use_result_type in c("tpr", "fpr")) {

  plot_df <- NULL
  
  # Check ALDEx2 error
  # Check scran error
  for(DE_method in c("DESeq2", "MAST", "NBGLM")) {
    
    # ----------------------------------------------------------------------------
    #   Call discrepancy by chosen method
    # ----------------------------------------------------------------------------
    
    # Get baseline differential abundance calls
    if(use_baseline == "threshold") {
      # "Differential" features will be those with mean fold change in abundance of 
      # >= 1.5 or <= 0.5
      oracle_calls <- calc_threshold_DA(ref_data)
    } else {
      oracle_calls <- NULL
    }
    
    rates <- calc_DE_discrepancy(ref_data,
                                 data,
                                 groups,
                                 method = DE_method,
                                 oracle_calls = oracle_calls)
    
    # ----------------------------------------------------------------------------
    #   Iterate TPR, FPR predictions
    # ----------------------------------------------------------------------------
    
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
      stop("Predictive model fit not found!")
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
    
    if(is.null(plot_df)) {
      pred_sim <- predict(fit_obj$result, newdata = fit_obj$test_features)
      
      plot_df <- data.frame(true = fit_obj$test_response,
                            predicted = pred_sim,
                            type = "simulated")
    }

    pred_real <- predict(fit_obj$result, newdata = features_df)
    
    plot_df <- rbind(plot_df,
                     data.frame(true = rates$tpr,
                                predicted = pred_real,
                                type = DE_method))
    
  }
  
  if(use_result_type == "tpr") {
    pl <- ggplot(plot_df, aes(x = true, y = predicted, fill = type))
  } else {
    pl <- ggplot(plot_df, aes(x = true, y = predicted, fill = type))
  }
  pl <- pl +
    geom_point(shape = 21, size = 3) +
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


