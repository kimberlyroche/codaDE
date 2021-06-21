source("path_fix.R")

library(codaDE)
library(tidyverse)
library(gridExtra)
library(randomForest)

# ------------------------------------------------------------------------------
#   Fit predictive model
# ------------------------------------------------------------------------------

# Fit predictive model
model_fn <- file.path("output", "RF_all.rds")
if(file.exists(model_fn)) {
  fit_obj <- readRDS(model_fn)
} else {
  fit_obj <- fit_predictive_model(model = "RF", do_predict = TRUE,
                                  DE_method = "all", plot_weights = FALSE)
  saveRDS(fit_obj, model_fn)
}

# ------------------------------------------------------------------------------
#   Parse and wrangle validation data
# ------------------------------------------------------------------------------

# dataset_name <- "Barlow"
dataset_name <- "Athanasiadou_ciona"
# dataset_name <- "Athanasiadou_yeast"

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

# Check scran pre-processing error here
use_method <- "DESeq2"

oracle_calls <- call_DA_NB(ref_data, groups)
oracle_calls <- oracle_calls$pval
res <- calc_DE_discrepancy(ref_data, data, groups, method = use_method,
                           oracle_calls = oracle_calls)

# res <- calc_DE_discrepancy(ref_data, data, groups, method = use_method, oracle_calls = oracle_calls)

cat(paste0("Baseline differential feature percentage: ",
           round(sum(oracle_calls < 0.05) / length(oracle_calls), 3)*100,
           "\n"))

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
features$CORRP <- 1 # Correlated features? (probably)
features$CORRP <- factor(features$CORRP, levels = c(0, 1))
features$PARTIAL <- 0 # No partial
features$PARTIAL <- factor(features$PARTIAL, levels = c(0, 1))
features$METHOD <- use_method
features$METHOD <- factor(features$METHOD, levels = c("ALDEx2", "DESeq2", "MAST", "NBGLM", "scran"))

# Convert to data.frame
features_df <- as.data.frame(features)

# Remove some features we're no longer using (too correlated with others)
features_df <- features_df %>%
  select(-c(FW_RA_PFC1_D, FW_CLR_MED_D, FW_CLR_SD_D, FW_CLR_PNEG_D))

plot_labels <- list(fpr = "specificity (1 - FPR)", tpr = "sensitivity (TPR)")

pl <- list()
for(use_result_type in c("tpr", "fpr")) {
  # Visualize this new data point in the mix
  test_response <- fit_obj$predictions[[use_result_type]]$true
  prediction <- fit_obj$predictions[[use_result_type]]$predicted
  p_labels <- fit_obj$predictions[[use_result_type]]$p_labels

  # Make predictions on this extra data point
  predict_ds <- predict(fit_obj$fitted_model[[use_result_type]], newdata = features_df)
  point_pred <- as.numeric(predict_ds)
  point_res <- as.numeric(res[[use_result_type]])
  if(use_result_type == "fpr") {
    # Convert to specificity
    point_res <- 1 - point_res
  }

  # Visualize TPR prediction
  plot_df <- data.frame(true = fit_obj$predictions[[use_result_type]]$true,
                        predicted = fit_obj$predictions[[use_result_type]]$predicted,
                        p = factor(fit_obj$predictions[[use_result_type]]$p_labels))
  levels(plot_df$p) <- c("low", "med", "high")
  pl[[use_result_type]] <- ggplot() +
    geom_point(data = plot_df,
               mapping = aes(x = true, y = predicted),
               shape = 21,
               size = 3,
               fill = "#bbbbbb") +
    geom_point(data = data.frame(true = point_pred,
                                 predicted = point_res),
               mapping = aes(x = true, y = predicted),
               shape = 21,
               size = 3,
               fill = "#ff0000") +
    xlim(c(0,1)) +
    ylim(c(0,1)) +
    labs(x = paste0("observed ", plot_labels[[use_result_type]]),
         y = paste0("predicted ", plot_labels[[use_result_type]]),
         fill = "feature no.")
}

grid.arrange(grobs = pl, ncol = 2)

# ------------------------------------------------------------------------------
#   Visualize features importance
# ------------------------------------------------------------------------------

use_result_type <- "tpr"

f_imp <- importance(fit_obj$fitted_model[[use_result_type]])

f_imp <- data.frame(feature = rownames(f_imp),
                   importance = unname(unlist(f_imp[,1])))
f_imp <- f_imp %>%
  filter(!(feature %in% c("METHOD", "PARTIAL", "CORRP")))

# ggplot(f_imp, aes(x = reorder(feature, importance), y = importance)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +
#   labs(x = "feature",
#        y = "importance ('mean decrease in node impurity')")

top_features <- f_imp %>%
  arrange(desc(importance)) %>%
  slice(1:20) %>%
  pull(feature)

train_features <- fit_obj$train_features[[use_result_type]]
train_features <- train_features[,colnames(train_features) %in% top_features]

plot_df <- pivot_longer(train_features, cols = everything(), names_to = "feature", values_to = "value")
plot_df2 <- pivot_longer(features_df[,colnames(features_df) %in% top_features],
                         cols = everything(), names_to = "feature", values_to = "value")
plot_df$feature <- factor(plot_df$feature, levels = top_features)
plot_df2$feature <- factor(plot_df2$feature, levels = top_features)

ggplot() +
  geom_boxplot(data = plot_df,
               mapping = aes(x = feature, y = value)) +
  geom_point(data = plot_df2,
             mapping = aes(x = feature, y = value), size = 3, shape = 21, fill = "red") +
  labs(x = "feature",
       y = "training range") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
  

