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
              help = "data set to use",
              metavar = "character"),
  make_option(c("--threshold"),
              type = "numeric",
              default = "1",
              help = "minimum mean abundance to threshold on",
              metavar = "numeric"),
  make_option(c("--modelsubdir"),
              type = "character",
              default = "regression",
              help = "predictive model sub directory",
              metavar = "character"),
  make_option(c("--permodel"),
              type = "logical",
              default = "TRUE",
              help = "use per-model predictive fits",
              metavar = "logical"),
  make_option(c("--selfbaseline"),
              type = "logical",
              default = "FALSE",
              help = "use self calls as reference instead of oracle",
              metavar = "logical"),
  make_option(c("--usetotals"),
              type = "logical",
              default = "FALSE",
              help = "use information about change in total abundances from absolute counts",
              metavar = "logical"),
  make_option(c("--userenormcounts"),
              type = "logical",
              default = "FALSE",
              help = "use features generated from renormalized counts",
              metavar = "logical"),
  make_option(c("--usecpm"),
              type = "logical",
              default = "FALSE",
              help = "scale relative abundances to counts per million",
              metavar = "logical")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

dataset_name <- opt$dataset
threshold <- opt$threshold
model_subdir <- opt$modelsubdir
per_model <- opt$permodel
use_self_baseline <- opt$selfbaseline
use_totals <- opt$usetotals
use_renorm_counts <- opt$userenormcounts
use_cpm <- opt$usecpm

testing <- FALSE
alpha <- 0.05
methods_list <- c("ALDEx2", "DESeq2", "scran")

if(threshold < 0) {
  stop(paste0("Invalid threshold: ", threshold, "!\n"))
}

# ------------------------------------------------------------------------------
#   Parse and wrangle validation data
# ------------------------------------------------------------------------------

abs_data <- do.call(paste0("parse_", dataset_name), list(absolute = TRUE))
rel_data <- do.call(paste0("parse_", dataset_name), list(absolute = FALSE))

print("A")

if(use_cpm) {
  # Convert to CPM
  rel_data$counts <- apply(rel_data$counts, 2, function(x) x/sum(x))
  rel_data$counts <- round(rel_data$counts*1e06)
}

if(testing & nrow(abs_data$counts) > 500) {
  k <- 500
  sample_idx <- sample(1:nrow(abs_data$counts), size = k, replace = FALSE)
  abs_data$counts <- abs_data$counts[sample_idx,]
  rel_data$counts <- rel_data$counts[sample_idx,]
}

print("B")

groups <- abs_data$groups
# Subsample if tons of cells/samples
downsample_limit <- 100 # was 100
set.seed(100)
pairs <- table(groups)
if(pairs[1] > downsample_limit) {
  A_sample <- sample(which(groups == names(pairs)[1]), size = downsample_limit)
  # A_sample <- which(groups == names(pairs)[1])[1:downsample_limit]
} else {
  A_sample <- which(groups == names(pairs)[1])
}
if(pairs[2] > downsample_limit) {
  B_sample <- sample(which(groups == names(pairs)[2]), size = downsample_limit)
  # B_sample <- which(groups == names(pairs)[2])[1:downsample_limit]
} else {
  B_sample <- which(groups == names(pairs)[2])
}
ref_data <- abs_data$counts[,c(A_sample, B_sample)]
data <- rel_data$counts[,c(A_sample, B_sample)]
groups <- groups[c(A_sample, B_sample)]

# Reorient as (samples x features)
ref_data <- t(ref_data)
data <- t(data)
groups <- factor(groups)
tax <- abs_data$tax

print("C")

# Convert to integers, just for DESeq2
ref_data <- apply(ref_data, c(1,2), as.integer)
data <- apply(data, c(1,2), as.integer)

print("D")

# Look at sparsity; filter out too-low features
retain_features <- colMeans(ref_data) >= threshold & colMeans(data) >= threshold
ref_data <- ref_data[,retain_features]
data <- data[,retain_features]
if(!is.null(tax)) {
  tax <- tax[retain_features]
}

print("E")

# ------------------------------------------------------------------------------
#   Call discrepancy by each method; save to file if doesn't already exist
#
#   This can be time-consuming!
# ------------------------------------------------------------------------------

# Save this estimate for later; we'll print it out with other statistics
percent_DE <- NULL
for(DE_method in methods_list) {
  # Pull saved calls on this data set x method if these exist
  save_fn <- file.path("output",
                       "real_data_calls",
                       "no_norm",
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
  if(is.null(percent_DE)) {
    oracle_calls <- p.adjust(readRDS(save_fn)$all_calls$oracle_calls, method = "BH") < 0.05
    percent_DE <- sum(oracle_calls) / length(oracle_calls)
  }
}

# ------------------------------------------------------------------------------
#   Wrangle data for prediction-making
# ------------------------------------------------------------------------------

# Get info we need from data to make a prediction
if(dataset_name == "VieiraSilva") {
  counts_A <- data[groups == "mHC",]
  counts_B <- data[groups == "CD",]
  counts_A_abs <- ref_data[groups == "mHC",]
  counts_B_abs <- ref_data[groups == "CD",]
}
if(dataset_name == "Barlow") {
  counts_A <- data[groups == "control",]
  counts_B <- data[groups == "keto",]
  counts_A_abs <- ref_data[groups == "control",]
  counts_B_abs <- ref_data[groups == "keto",]
}
if(dataset_name == "Song") {
  counts_A <- data[groups == "lung",]
  counts_B <- data[groups == "brain",]
  counts_A_abs <- ref_data[groups == "lung",]
  counts_B_abs <- ref_data[groups == "brain",]
}
if(dataset_name == "Monaco") {
  counts_A <- data[groups == "CD4_naive",]
  counts_B <- data[groups == "PBMC",]
  counts_A_abs <- ref_data[groups == "CD4_naive",]
  counts_B_abs <- ref_data[groups == "PBMC",]
}
if(dataset_name == "Hagai") {
  counts_A <- data[groups == "unstimulated",]
  counts_B <- data[groups == "pIC4",]
  counts_A_abs <- ref_data[groups == "unstimulated",]
  counts_B_abs <- ref_data[groups == "pIC4",]
}
if(dataset_name == "Owens") {
  counts_A <- data[groups == "early_series",]
  counts_B <- data[groups == "late_series",]
  counts_A_abs <- ref_data[groups == "early_series",]
  counts_B_abs <- ref_data[groups == "late_series",]
}
if(dataset_name == "Klein") {
  counts_A <- data[groups == "unstimulated",]
  counts_B <- data[groups == "LIF-2hr",]
  counts_A_abs <- ref_data[groups == "unstimulated",]
  counts_B_abs <- ref_data[groups == "LIF-2hr",]
}
if(dataset_name == "Yu") {
  counts_A <- data[groups == "Brn",]
  counts_B <- data[groups == "Lvr",]
  counts_A_abs <- ref_data[groups == "Brn",]
  counts_B_abs <- ref_data[groups == "Lvr",]
}
if(dataset_name == "Muraro") {
  counts_A <- data[groups == "alpha",]
  counts_B <- data[groups == "beta",]
  counts_A_abs <- ref_data[groups == "alpha",]
  counts_B_abs <- ref_data[groups == "beta",]
}
if(dataset_name == "Hashimshony") {
  counts_A <- data[groups == "0",]
  counts_B <- data[groups == "1",]
  counts_A_abs <- ref_data[groups == "0",]
  counts_B_abs <- ref_data[groups == "1",]
}
if(dataset_name == "Kimmerling") {
  counts_A <- data[groups == "low_mass",]
  counts_B <- data[groups == "high_mass",]
  counts_A_abs <- ref_data[groups == "low_mass",]
  counts_B_abs <- ref_data[groups == "high_mass",]
}
if(dataset_name == "Gruen") {
  counts_A <- data[groups == "A",]
  counts_B <- data[groups == "B",]
  counts_A_abs <- ref_data[groups == "A",]
  counts_B_abs <- ref_data[groups == "B",]
}

print("G")

save_fn <- file.path("output",
                     paste0("filtered_data_",dataset_name,"_threshold",threshold,".rds"))
if(!file.exists(save_fn)) {
  saveRDS(list(absolute = ref_data,
               relative = data,
               groups = groups), save_fn)
}

# Report stats on this data set
cat(paste0("Number of features (after filtering): ", ncol(counts_A), "\n"))
cat(paste0("Number of samples: ", length(groups), "\n"))
cat(paste0("Samples per condition (A, B): ", sum(groups == levels(groups)[1]), ", ",
           sum(groups == levels(groups)[2]), "\n"))
cat(paste0("Percent zeros: ", round(sum(data == 0)/(nrow(data)*ncol(data)), 3)*100, "%\n"))
sA <- mean(rowSums(counts_A_abs))
sB <- mean(rowSums(counts_B_abs))
fc <- max(sA, sB) / min(sA, sB)
cat(paste0("Approx. fold change between conditions: ", round(fc, 1), "\n"))
cat(paste0("Approx. percent differential features: ", round(percent_DE, 2)*100, "%\n"))

print("H")

# This takes 2-3 min. to run on 15K features
features <- as.data.frame(characterize_dataset(counts_A, counts_B))

print("I")

save_fn <- file.path("output",
                     paste0("filtered_features_",dataset_name,"_threshold",threshold,".rds"))
if(!file.exists(save_fn)) {
  saveRDS(features, save_fn)
}

features <- features %>%
  select(-c(FW_RA_PFC1_D, FW_CLR_MED_D, FW_CLR_SD_D, FW_CLR_PNEG_D))

# Add a few more
features$P <- ncol(counts_A)

if(use_totals) {
  m1 <- mean(rowSums(counts_A_abs))
  m2 <- mean(rowSums(counts_B_abs))
  m2 / m1
  features <- cbind(features, FC_ABSOLUTE = m2 / m1)
}

if(use_renorm_counts) {
  scaled_DESeq2 <- scaled_counts_DESeq2(data, groups, pseudocount = 0.5)
  features_DESeq2 <- as.data.frame(characterize_dataset(scaled_DESeq2[groups == levels(groups)[1],],
                                                        scaled_DESeq2[groups == levels(groups)[2],])) %>%
    select(-c(FW_RA_PFC1_D, FW_CLR_MED_D, FW_CLR_SD_D, FW_CLR_PNEG_D))
  
  ignore_columns <- c("P", "FC_ABSOLUTE")
  rename_flag <- !(colnames(features) %in% ignore_columns)
  colnames(features)[rename_flag] <- paste0(colnames(features)[rename_flag], "_rel")
  rename_flag <- !(colnames(features_DESeq2) %in% ignore_columns)
  colnames(features_DESeq2)[rename_flag] <- paste0(colnames(features_DESeq2)[rename_flag], "_scaled")
  
  features <- cbind(features, features_DESeq2)
}

# ------------------------------------------------------------------------------
#   Make predictions on simulations and real data and visualize these together
# ------------------------------------------------------------------------------

print("J")

save_df <- NULL

for(use_result_type in c("TPR", "FPR")) {

  print("K")

  # Load predictive model
  model_list <- c("all")
  if(per_model) {
    model_list <- methods_list
  }
  for(method_type in model_list) {
    if(method_type == "all") {
      model_fn <- file.path("output",
                            "predictive_fits",
                            ifelse(use_self_baseline, "self", "oracle"),
                            model_subdir,
                            paste0(use_result_type, ".rds"))
    } else {
      model_fn <- file.path("output",
                            "predictive_fits",
                            ifelse(use_self_baseline, "self", "oracle"),
                            model_subdir,
                            paste0(use_result_type,
                                   "_",
                                   method_type,
                                   ".rds"))
    }
    if(!file.exists(model_fn)) {
      stop(paste0("Predictive model fit not found: ", model_fn, "!\n"))
    }
    fit_obj <- readRDS(model_fn)
  
    for(DE_method in methods_list) {
      
      # Skip this method if we're using model-specific predictions and this
      # iteration doesn't match the outer loop
      if(per_model & DE_method != method_type) {
        next;
      }
      
      # Skip this method if we didn't train on it!
      if(!per_model &
         !is.null(fit_obj$result$forest$xlevels$METHOD) &
         !(DE_method %in% fit_obj$result$forest$xlevels$METHOD)) {
        next;
      }
  
      cat(paste0("Evaluating ", use_result_type, " x ", DE_method, "\n"))
  
      # --------------------------------------------------------------------------
      #   Pull calls for chosen method
      # --------------------------------------------------------------------------
      
      save_fn <- file.path("output",
                           "real_data_calls",
                           "no_norm",
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
  
      # --------------------------------------------------------------------------
      #   Make predictions
      # --------------------------------------------------------------------------
  
      if(!per_model) {
        features$METHOD <- DE_method
        features$METHOD <- factor(features$METHOD,
                                  levels = fit_obj$result$forest$xlevels$METHOD)
      }

      pred_real <- predict(fit_obj$result, newdata = features, predict.all = TRUE)

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

save_path <- list("output",
                  "predictive_fits",
                  ifelse(use_self_baseline, "self", "oracle"),
                  model_subdir,
                  "validation_results",
                  "no_norm")
for(i in 1:length(save_path)) {
  fp <- do.call(file.path, save_path[1:i])
  if(!dir.exists(fp)) {
    dir.create(fp)
  }
}
save_fn <- file.path(fp,
                     paste0("results_",
                     dataset_name,
                     "_threshold",
                     threshold,
                     ".tsv"))
cat(paste0("Saving output to: ", save_fn, "\n"))
write.table(save_df, save_fn, sep = "\t", quote = FALSE, row.names = FALSE)
