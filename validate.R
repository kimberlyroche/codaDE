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
  make_option(c("--baseline"),
              type = "character",
              default = "self",
              help = "calls to use as a reference: self, oracle",
              metavar = "character"),
  make_option(c("--threshold"),
              type = "numeric",
              default = "5",
              help = "minimum mean abundance to threshold on",
              metavar = "numeric"),
  make_option(c("--model_folder"),
              type = "character",
              default = "self_nopartial",
              help = "folder containing fitted models from which to predict",
              metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

dataset_name <- opt$dataset
use_baseline <- opt$baseline
threshold <- opt$threshold
model_dir <- opt$model_folder
testing <- FALSE

methods_list <- c("ALDEx2", "DESeq2", "scran")

if(!(dataset_name %in% c("VieiraSilva", "Barlow", "Song",
                         "Monaco", "Hagai", "Owens", "Klein", "Yu"))) {
  stop(paste0("Invalid data set: ", dataset_name, "!\n"))
}

if(!(use_baseline %in% c("self", "oracle"))) {
  stop(paste0("Invalid DA baseline: ", use_baseline, "!\n"))
}

if(threshold < 0) {
  stop(paste0("Invalid threshold: ", threshold, "!\n"))
}

if(!file.exists(file.path("output", "predictive_fits", model_dir))) {
  stop(paste0("Model folder does not exist: ", model_dir, "\n"))
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
# features_simulated <- pull_features(use_baseline = use_baseline,
#                                     feature_list = c("METHOD"))
# features_simulated <- features_simulated %>%
#     filter(METHOD %in% methods_list)
# features_simulated$METHOD <- factor(features_simulated$METHOD)

save_df <- NULL

for(use_result_type in c("TPR", "FPR")) {

  # Load predictive model
  model_fn <- file.path("output",
                        "predictive_fits",
                        model_dir,
                        paste0(use_baseline,
                               "_",
                               use_result_type,
                               ".rds"))
  if(!file.exists(model_fn)) {
    stop(paste0("Predictive model fit not found: ", model_fn, "!\n"))
  }
  fit_obj <- readRDS(model_fn)

  for(DE_method in methods_list) {
    
    if(!(DE_method %in% fit_obj$result$forest$xlevels$METHOD)) {
      next;
    }

    cat(paste0("Evaluating ", use_result_type, " x ", DE_method, "\n"))

    # --------------------------------------------------------------------------
    #   Pull calls for chosen method
    # --------------------------------------------------------------------------
    
    save_fn <- file.path("output",
                         "real_data_calls",
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
      stop(paste0("Missing file ", save_fn, "!"))
    }
    calls_obj <- readRDS(save_fn)
    all_calls <- calls_obj$all_calls
    rates <- calls_obj$rates
    rm(calls_obj)

    # --------------------------------------------------------------------------
    #   Make predictions
    # --------------------------------------------------------------------------

    save_fn <- file.path("output",
                         "predictive_fits",
                         model_dir,
                         "real_data_predictions",
                         paste0("predictions_",
                                use_baseline,
                                "_",
                                DE_method,
                                "_",
                                dataset_name,
                                "_threshold",
                                threshold,
                                ".rds"))
    if(file.exists(save_fn)) {
      pred_real <- readRDS(save_fn)
    } else {
      features$METHOD <- DE_method
      features$METHOD <- factor(features$METHOD,
                                levels = fit_obj$result$forest$xlevels$METHOD)
  
      pred_real <- predict(fit_obj$result, newdata = features, predict.all = TRUE)
      saveRDS(pred_real, save_fn)
    }
    
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

save_fn <- file.path("output",
                     "predictive_fits",
                     model_dir,
                     paste0("results_",
                            dataset_name,
                            "_threshold",
                            threshold,
                            ".tsv"))
write.table(save_df, save_fn, sep = "\t", quote = FALSE, row.names = FALSE)

# LEFT OFF HERE
# Need to add an FC_ABSOLUTE estimate and/or PERCENT_DE estimate for the real data feature set!
# Just base this on model_dir
