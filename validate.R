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
              help = "data set to use: VieiraSilva, Barlow, Song, Monaco, Muraro, Hashimshony, Kimmerling",
              metavar = "character"),
  make_option(c("--threshold"),
              type = "numeric",
              default = "5",
              help = "minimum mean abundance to threshold on",
              metavar = "numeric"),
  make_option(c("--norm"),
              type = "logical",
              default = "FALSE",
              help = "normalize observed total abundances",
              metavar = "logical"),
  make_option(c("--classify"),
              type = "logical",
              default = "FALSE",
              help = "classification (vs. regression) flag",
              metavar = "logical"),
  make_option(c("--permodel"),
              type = "logical",
              default = "TRUE",
              help = "use per-model predictive fits",
              metavar = "logical"),
  make_option(c("--selfbaseline"),
              type = "logical",
              default = "TRUE",
              help = "use self calls as reference instead of oracle",
              metavar = "logical"),
  make_option(c("--partials"),
              type = "logical",
              default = "FALSE",
              help = "flag indicating whether or not to use simulations with partially informative total abundances",
              metavar = "logical")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

dataset_name <- opt$dataset
threshold <- opt$threshold
do_norm <- opt$norm
do_classify <- opt$classify
per_model <- opt$permodel

use_self_baseline <- opt$selfbaseline
use_partials <- opt$partials

testing <- FALSE

methods_list <- c("ALDEx2", "DESeq2", "scran")

if(!(dataset_name %in% c("VieiraSilva", "Barlow", "Song",
                         "Monaco", "Hagai", "Owens", "Klein", "Yu"))) {
  stop(paste0("Invalid data set: ", dataset_name, "!\n"))
}

if(threshold < 0) {
  stop(paste0("Invalid threshold: ", threshold, "!\n"))
}

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
downsample_limit <- 100 # was 100
set.seed(100)
pairs <- table(groups)
if(pairs[1] > downsample_limit) {
  A_sample <- sample(which(groups == names(pairs)[1]), size = downsample_limit)
} else {
  A_sample <- which(groups == names(pairs)[1])
}
if(pairs[2] > downsample_limit) {
  B_sample <- sample(which(groups == names(pairs)[2]), size = downsample_limit)
} else {
  B_sample <- which(groups == names(pairs)[2])
}
ref_data <- ref_data[c(A_sample, B_sample),]
data <- data[c(A_sample, B_sample),]
groups <- groups[c(A_sample, B_sample)]

# Convert to integers, just for DESeq2
ref_data <- apply(ref_data, c(1,2), as.integer)
data <- apply(data, c(1,2), as.integer)

# New 2021/09/09: normalize relative abundances
# Each column should sum to roughly the same abundance
if(do_norm) {
  new_total_scale <- median(rowSums(data))
  data <- t(apply(data, 1, function(x) round((x/sum(x))*new_total_scale)))
}

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

for(DE_method in methods_list) {
  # Pull saved calls on this data set x method if these exist
  save_fn <- file.path("output",
                       "real_data_calls",
                       ifelse(do_norm, "norm", "no_norm"),
                       paste0("calls_",
                              ifelse(use_self_baseline, "self", "oracle"),
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
    if(use_self_baseline) {
      oracle_calls <- NULL
    } else {
      oracle_calls <- call_DA_NB(ref_data, groups)$pval
    }
    all_calls <- DA_wrapper(ref_data, data, groups, DE_method, oracle_calls)
    rates <- calc_DA_discrepancy(all_calls$calls, all_calls$oracle_calls)
    save_obj <- list(all_calls = all_calls, rates = rates)
    saveRDS(save_obj, save_fn)
  }
}

# ------------------------------------------------------------------------------
#   Wrangle data for prediction-making
# ------------------------------------------------------------------------------

# Get info we need from data to make a prediction
if(dataset_name == "VieiraSilva") {
  counts_A <- data[groups == "mHC",]
  counts_B <- data[groups == "CD",]
}
if(dataset_name == "Barlow") {
  counts_A <- data[groups == "control",]
  counts_B <- data[groups == "keto",]
}
if(dataset_name == "Song") {
  counts_A <- data[groups == "lung",]
  counts_B <- data[groups == "brain",]
}
if(dataset_name == "Monaco") {
  counts_A <- data[groups == "CD4_naive",]
  counts_B <- data[groups == "PBMC",]
}
if(dataset_name == "Hagai") {
  counts_A <- data[groups == "unstimulated",]
  counts_B <- data[groups == "pIC4",]
}
if(dataset_name == "Owens") {
  counts_A <- data[groups == "early_series",]
  counts_B <- data[groups == "late_series",]
}
if(dataset_name == "Klein") {
  counts_A <- data[groups == "unstimulated",]
  counts_B <- data[groups == "LIF-2hr",]
}
if(dataset_name == "Yu") {
  counts_A <- data[groups == "Brn",]
  counts_B <- data[groups == "Lvr",]
}

# This takes 2-3 min. to run on 15K features
features <- as.data.frame(characterize_dataset(counts_A, counts_B))

# Add a few more
features$P <- ncol(counts_A)

features <- features %>%
  select(c(-FW_RA_PFC1_D, FW_CLR_MED_D, FW_CLR_SD_D, FW_CLR_PNEG_D))

# ------------------------------------------------------------------------------
#   Make predictions on simulations and real data and visualize these together
# ------------------------------------------------------------------------------

# It's crucial we match the factor levels the model was trained on
# Features will return all levels for METHOD (including methods we evaluated
# on when testing but aren't reporting). Make sure these 
# features_simulated <- pull_features(use_baseline = ifelse(use_self_baseline, "self", "oracle"),
#                                     feature_list = c("METHOD"))
# features_simulated <- features_simulated %>%
#     filter(METHOD %in% methods_list)
# features_simulated$METHOD <- factor(features_simulated$METHOD)

save_df <- NULL

for(use_result_type in c("TPR", "FPR")) {

  # Load predictive model
  model_list <- c("all")
  if(per_model) {
    model_list <- methods_list
  }
  for(model_type in model_list) {
    if(model_type == "all") {
      model_fn <- file.path("output",
                            "predictive_fits",
                            paste0(ifelse(use_self_baseline, "self", "oracle"),
                                   "_",
                                   ifelse(use_partials, "partial", "nopartial")),
                            ifelse(do_classify, "classification", "regression"),
                            paste0(use_result_type, ".rds"))
    } else {
      model_fn <- file.path("output",
                            "predictive_fits",
                            paste0(ifelse(use_self_baseline, "self", "oracle"),
                                   "_",
                                   ifelse(use_partials, "partial", "nopartial")),
                            ifelse(do_classify, "classification", "regression"),
                            paste0(use_result_type,
                                   "_",
                                   model_type,
                                   ".rds"))
    }
    if(!file.exists(model_fn)) {
      stop(paste0("Predictive model fit not found: ", model_fn, "!\n"))
    }
    fit_obj <- readRDS(model_fn)
  
    for(DE_method in methods_list) {
      
      # Skip this method if we're using model-specific predictions and this
      # iteration doesn't match the outer loop
      if(per_model & DE_method != model_type) {
        next;
      }
      
      # Skip this method if we didn't train on it!
      if(!per_model & !(DE_method %in% fit_obj$result$forest$xlevels$METHOD)) {
        next;
      }
  
      cat(paste0("Evaluating ", use_result_type, " x ", DE_method, "\n"))
  
      # --------------------------------------------------------------------------
      #   Pull calls for chosen method
      # --------------------------------------------------------------------------
      
      save_fn <- file.path("output",
                           "real_data_calls",
                           ifelse(do_norm, "norm", "no_norm"),
                           paste0("calls_",
                                  ifelse(use_self_baseline, "self", "oracle"),
                                  "_",
                                  DE_method,
                                  "_",
                                  dataset_name,
                                  "_threshold",
                                  threshold,
                                  ".rds"))
      if(!file.exists(save_fn)) {
        stop(paste0("Missing file ", save_fn, "!"))
      }
      calls_obj <- readRDS(save_fn)
      all_calls <- calls_obj$all_calls
      rates <- calls_obj$rates
      rm(calls_obj)

      if(use_result_type == "TPR") {
        true_result <- ifelse(rates$TPR < 0.95, 1, 0)
      } else {
        true_result <- ifelse(1 - rates$FPR < 0.95, 1, 0)
      }
      true_result <- factor(true_result, levels = c(0, 1))
  
      # --------------------------------------------------------------------------
      #   Make predictions
      # --------------------------------------------------------------------------
  
      if(!per_model) {
        features$METHOD <- DE_method
        features$METHOD <- factor(features$METHOD,
                                  levels = fit_obj$result$forest$xlevels$METHOD)
      }

      pred_real <- predict(fit_obj$result, newdata = features, predict.all = TRUE)

      if(do_classify) {
        # Just use the aggregate prediction
        save_df <- rbind(save_df,
                         data.frame(dataset = dataset_name,
                                    result_type = "true",
                                    DE_method = DE_method,
                                    threshold = threshold,
                                    score_type = use_result_type,
                                    point = true_result))
        save_df <- rbind(save_df,
                         data.frame(dataset = dataset_name,
                                    result_type = "predicted",
                                    DE_method = DE_method,
                                    threshold = threshold,
                                    score_type = use_result_type,
                                    point = pred_real$aggregate))
      } else {        
        save_df <- rbind(save_df,
                         data.frame(dataset = dataset_name,
                                    result_type = "true",
                                    DE_method = DE_method,
                                    threshold = threshold,
                                    score_type = use_result_type,
                                    lower90 = NA,
                                    lower50 = NA,
                                    upper50 = NA,
                                    upper90 = NA,
                                    point = ifelse(use_result_type == "TPR",
                                                   rates$TPR,
                                                   1 - rates$FPR)))
        save_df <- rbind(save_df,
                         data.frame(dataset = dataset_name,
                                    result_type = "predicted",
                                    DE_method = DE_method,
                                    threshold = threshold,
                                    score_type = use_result_type,
                                    lower90 = unname(quantile(pred_real$individual[1,], probs = c(0.05))),
                                    lower50 = unname(quantile(pred_real$individual[1,], probs = c(0.25))),
                                    upper50 = unname(quantile(pred_real$individual[1,], probs = c(0.75))),
                                    upper90 = unname(quantile(pred_real$individual[1,], probs = c(0.95))),
                                    point = pred_real$aggregate))
      }
    }
  }
}

save_fn <- file.path("output",
                     "predictive_fits",
                     paste0(ifelse(use_self_baseline, "self", "oracle"),
                                   "_",
                                   ifelse(use_partials, "partial", "nopartial")),
                     ifelse(do_classify, "classification", "regression"),
                     "validation_results",
                     ifelse(do_norm, "norm", "no_norm"),
                     paste0("results_",
                            dataset_name,
                            "_threshold",
                            threshold,
                            ".tsv"))
cat(paste0("Saving output to: ", save_fn, "\n"))
write.table(save_df, save_fn, sep = "\t", quote = FALSE, row.names = FALSE)
