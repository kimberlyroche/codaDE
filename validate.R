source("path_fix.R")

library(codaDE)
library(tidyverse)
library(optparse)

# ------------------------------------------------------------------------------
#   All this (from `evaluate_methods.R`) needs to go into permanent functions!
# ------------------------------------------------------------------------------

DE_by_DESeq2 <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_DESeq2(ref_data, groups)
  }
  calls <- call_DA_DESeq2(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

calc_DE_discrepancy <- function(ref_data, data, groups, method = "NBGLM",
                                oracle_calls = NULL) {
  if(method == "NBGLM") {
    DE_calls <- DE_by_NB(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
  if(method == "DESeq2") {
    DE_calls <- DE_by_DESeq2(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
  if(method == "MAST") {
    DE_calls <- DE_by_MAST(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
  if(method == "ALDEx2") {
    DE_calls <- DE_by_ALDEx2(ref_data, data, groups,
                             oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
  if(method == "scran") {
    DE_calls <- DE_by_scran(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
  
  de <- calls < 0.05
  sim_de <- oracle_calls < 0.05
  
  TP <- sum(de & sim_de)
  FP <- sum(de & !sim_de)
  
  TN <- sum(!de & !sim_de)
  FN <- sum(!de & sim_de)
  
  return(list(tpr = TP/(TP+FN), fpr = FP/(FP+TN), oracle_calls = oracle_calls))
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

abs_data <- parse_Barlow(absolute = TRUE)
rel_data <- parse_Barlow(absolute = FALSE)

# Get discrepancy
ref_data <- t(abs_data$counts)
data <- t(rel_data$counts)
groups <- abs_data$groups

# Convert to integers, just for DESeq2
ref_data <- apply(ref_data, c(1,2), as.integer)
data <- apply(data, c(1,2), as.integer)

res <- calc_DE_discrepancy(ref_data, data, groups, method = "DESeq2", oracle_calls = NULL)

# Get info we need from data to make a prediction
counts_A <- rel_data$counts[,rel_data$groups == "control"]
counts_B <- rel_data$counts[,rel_data$groups == "keto"]

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

features <- characterize_dataset(t(counts_A), t(counts_B))

# Add a few more
features$P <- nrow(counts_A)
features$CORRP <- 1 # Correlated features? (probably)
features$CORRP <- factor(features$CORRP, levels = c(0, 1))
features$PARTIAL <- 0 # No partial
features$PARTIAL <- factor(features$PARTIAL, levels = c(0, 1))
features$METHOD <- "DESeq2"
features$METHOD <- factor(features$METHOD, levels = c("ALDEx2", "DESeq2", "MAST", "NBGLM", "scran"))

# Convert to data.frame
features_df <- as.data.frame(features)

# Remove some features we're no longer using (too correlated with others)
features_df <- features_df %>%
  select(-c(FW_RA_PFC1_D, FW_CLR_MED_D, FW_CLR_SD_D, FW_CLR_PNEG_D))

# Fit predictive model
fit_obj <- fit_predictive_model(model = "RF", do_predict = TRUE,
                                DE_method = "all", plot_weights = FALSE)

# Visualize this new data point in the mix
test_response <- fit_obj$predictions[[use_result_type]]$true
prediction <- fit_obj$predictions[[use_result_type]]$predicted
p_labels <- fit_obj$predictions[[use_result_type]]$p_labels

for(use_result_type in c("tpr", "fpr")) {
  # Make predictions on this extra data point
  predict_ds <- predict(fit_obj$fitted_model[[use_result_type]], newdata = features_df)
  if(use_result_type == "tpr") {
    point_pred <- as.numeric(predict_ds)
  } else {
    point_pred <- 1 - as.numeric(predict_ds)
  }
  
  # Visualize TPR prediction
  plot_df <- data.frame(true = fit_obj$predictions[[use_result_type]]$true,
                        predicted = fit_obj$predictions[[use_result_type]]$predicted,
                        p = factor(fit_obj$predictions[[use_result_type]]$p_labels))
  levels(plot_df$p) <- c("low", "med", "high")
  pl <- ggplot() +
    geom_point(data = plot_df,
               mapping = aes(x = true, y = predicted),
               shape = 21,
               size = 3,
               fill = "#bbbbbb") +
    geom_point(data = data.frame(true = point_pred,
                                 predicted = res[[use_result_type]]),
               mapping = aes(x = true, y = predicted),
               shape = 21,
               size = 3,
               fill = "#ff0000") +
    xlim(c(0,1)) +
    ylim(c(0,1)) +
    labs(x = paste0("observed ", plot_labels[[use_result_type]]),
         y = paste0("predicted ", plot_labels[[use_result_type]]),
         fill = "feature no.")
  show(pl)
}

