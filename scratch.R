source("path_fix.R")

library(tidyverse)
library(codaDE)
library(RSQLite)
library(gridExtra)
library(RColorBrewer)
library(randomForest)
library(entropy)

source("ggplot_fix.R")

# -----------------------------------------------------------------------------
#   Functions
# -----------------------------------------------------------------------------

get_majority_sign <- function(counts_A, counts_B) {
  m1 <- mean(rowSums(counts_A)) # sample-wise
  m2 <- mean(rowSums(counts_B))
  sign(m2 - m1)
}

get_featurewise_change <- function(counts_A, counts_B) {
  m1 <- colMeans(counts_A) # feature-wise
  m2 <- colMeans(counts_B)
  m2 - m1
}

get_prop_sign_agreement <- function(counts_A, counts_B) {
  # How many features change in the same direction as overall change?
  majority_sign <- get_majority_sign(counts_A, counts_B)
  change <- get_featurewise_change(counts_A, counts_B)
  sum(sign(change) == majority_sign) / length(change)
}

get_prop_top10 <- function(counts_A, counts_B) {
  # The top 10% of features account for what % of the majority change?
  majority_sign <- get_majority_sign(counts_A, counts_B)
  change <- get_featurewise_change(counts_A, counts_B)
  majority_change <- change[sign(change) == majority_sign]
  ranks <- order(majority_change, decreasing = TRUE)
  n <- length(ranks)
  n10 <- round(n/10)
  sum(majority_change[ranks[1:n10]]) / sum(majority_change)
}

get_max_diff <- function(counts_A, counts_B) {
  # What's the max differential for a single feature as a fraction of majority change?
  majority_sign <- get_majority_sign(counts_A, counts_B)
  change <- get_featurewise_change(counts_A, counts_B)
  majority_change <- change[sign(change) == majority_sign]
  max(majority_change) / sum(majority_change)
}

get_max_diff_relab <- function(counts_A, counts_B) {
  # What's the max differential for a single feature as a fraction of majority change?
  m1 <- colMeans(counts_A)
  m1 <- m1 / sum(m1)
  m2 <- colMeans(counts_B)
  m2 <- m2 / sum(m2)
  max(abs(m2 - m1))
}

get_prop_stable <- function(counts_A, counts_B) {
  # Proportion "stable" features (i.e. non-negligibly abundant with less than
  # 1.5-fold change)
  m1 <- colMeans(counts_A) # feature-wise
  m2 <- colMeans(counts_B)
  include_features <- m2 >= 1 & m1 >= 1
  feature_fc <- m2[include_features] / m1[include_features]
  feature_fc[feature_fc < 1] <- 1 / feature_fc[feature_fc < 1]
  sum(feature_fc < 1.5) / length(feature_fc)
}

get_entropy <- function(counts_A, counts_B, bins = 20) {
  # Entropy of change
  change <- get_featurewise_change(counts_A, counts_B)
  entropy(table(cut(change, breaks = bins))) / entropy(rep(1, bins))
}

get_samples_agree <- function(counts_A, counts_B) {
  # How many sample pairs agree on the direction of change?
  m1 <- rowSums(counts_A)
  m2 <- rowSums(counts_B)
  diffs <- c()
  for(i in 1:length(m1)) {
    for(j in 1:length(m2)) {
      diffs <- c(diffs, m2[j] - m1[j])
    }
  }
  pos_n <- sum(sign(diffs) > 0)
  neg_n <- sum(sign(diffs) < 0)
  sum(max(pos_n, neg_n)) / length(diffs)
}

get_zeros <- function(counts_A, counts_B) {
  # What's the percent zeros?
  max(sum(counts_A == 0)/(nrow(counts_A)*ncol(counts_A)),
      sum(counts_B == 0)/(nrow(counts_A)*ncol(counts_A)))
}

get_library_variation <- function(counts_A, counts_B) {
  # What percent of the largest library size is the smallest library size?
  # I.e. how much variation in library size is there?
  m1 <- rowSums(counts_A)
  m2 <- rowSums(counts_B)
  min(min(m1) / max(m1), min(m2) / max(m2))
}

get_depth_fc <- function(absolute_counts, relative_counts) {
  absolute_total_scale <- median(rowSums(absolute_counts))
  relative_total_scale <- median(rowSums(relative_counts))
  relative_total_scale / absolute_total_scale
}

# uuid <- "b20e0f28-4f58-4338-98f1-b174f639d552" # bad
# uuid <- "4cdb8ed4-db95-41b3-8aa1-d0f2d71fb5ce" # good
# temp <- readRDS(paste0("output/datasets/", uuid, ".rds"))
# absolute_counts <- temp$simulation$abundances
# relative_counts <- temp$simulation$observed_counts1
# 
# plot_stacked_bars(absolute_counts)
# plot_stacked_bars(relative_counts)
# 
# get_depth_fc(absolute_counts, relative_counts)
# 
# plot(c(rowSums(absolute_counts), rowSums(relative_counts)))
# 
# conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
# res <- dbGetQuery(conn, paste0("SELECT ",
#                                "datasets.BASELINE_CALLS AS ORACLE_BASELINE, ",
#                                "CALLS, ",
#                                "results.BASELINE_CALLS AS SELF_BASELINE ",
#                                "FROM results LEFT JOIN datasets ON ",
#                                "results.UUID=datasets.UUID ",
#                                "WHERE results.UUID='", uuid, "' AND ",
#                                "METHOD='DESeq2' AND PARTIAL_INFO=0 AND BASELINE_TYPE='self'"))
# dbDisconnect(conn)
# 
# oracle_abs <- as.numeric(strsplit(res$ORACLE_BASELINE, ";")[[1]])
# self_abs <- as.numeric(strsplit(res$SELF_BASELINE, ";")[[1]])
# self_rel <- as.numeric(strsplit(res$CALLS, ";")[[1]])
# 
# oracle_abs <- p.adjust(oracle_abs, method = "BH")
# self_abs <- p.adjust(self_abs, method = "BH")
# self_rel <- p.adjust(self_rel, method = "BH")
# 
# sum(oracle_abs < 0.05)
# sum(self_abs < 0.05)
# sum(self_rel < 0.05)
# 
# TP <- sum(self_abs < 0.05 & self_rel < 0.05)
# TN <- sum(self_abs >= 0.05 & self_rel >= 0.05)
# FP <- sum(self_abs >= 0.05 & self_rel < 0.05)
# FN <- sum(self_abs < 0.05 & self_rel >= 0.05)
# 
# FP
# TN
# FP / (FP + TN)

# -----------------------------------------------------------------------------
#   Parse real data
# -----------------------------------------------------------------------------

datasets <- c("VieiraSilva", "Barlow", "Song", "Monaco", "Hagai", "Owens", "Klein", "Yu")
thresholds <- c(1, 1, 1, 2, 3, 1, 1, 1);
fprs <- 1 - c(0.76, 0.95, 0.95, 0.93, 0.98, 0.80, 0.90, 0.90); # DESeq2
# fprs <- 1 - c(1, 0.80, 0.98, 1, 1, 0.74, 0.80, 0.78); # ALDEx2
# fprs <- 1 - c(0.93, 0.97, 0.99, 0.97, 0.99, 0.82, 0.80, 0.93); # scran

save_fn <- "bad_good_realdatasets.rds"
if(file.exists(save_fn)) {
  real_datasets <- readRDS(save_fn)
} else {
  real_datasets <- list()

  for(i in 1:length(datasets)) {
    dataset_name <- datasets[i]
    threshold <- thresholds[i]
    cat(paste0("Parsing ", dataset_name, "\n"))

    abs_data <- do.call(paste0("parse_", dataset_name), list(absolute = TRUE))
    rel_data <- do.call(paste0("parse_", dataset_name), list(absolute = FALSE))

    # Reorient as (samples x features)
    ref_data <- t(abs_data$counts)
    data <- t(rel_data$counts)
    groups <- abs_data$groups
    groups <- factor(groups)

    # Subsample if tons of cells/samples
    downsample_limit <- 100
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
    # cat(paste0("Percent zeros: ", round(sum(data == 0)/(nrow(data)*ncol(data)), 3)*100, "%\n"))

    # Get info we need from data to make a prediction
    if(dataset_name == "VieiraSilva") {
      label1 <- "mHC"
      label2 <- "CD"
    } else if(dataset_name == "Barlow") {
      label1 <- "control"
      label2 <- "keto"
    } else if(dataset_name == "Song") {
      label1 <- "lung"
      label2 <- "brain"
    } else if(dataset_name == "Monaco") {
      label1 <- "CD4_naive"
      label2 <- "PBMC"
    } else if(dataset_name == "Hagai") {
      label1 <- "unstimulated"
      label2 <- "pIC4"
    } else if(dataset_name == "Owens") {
      label1 <- "early_series"
      label2 <- "late_series"
    } else if(dataset_name == "Klein") {
      label1 <- "unstimulated"
      label2 <- "LIF-2hr"
    } else if(dataset_name == "Yu") {
      label1 <- "Brn"
      label2 <- "Lvr"
    } else {
      stop("Unknown data set name!")
    }
    real_datasets[[dataset_name]] <- list(absolute = list(ref_data[groups == label1,],
                                                          ref_data[groups == label2,]),
                                          relative = list(data[groups == label1,],
                                                          data[groups == label2,]))
  }
  saveRDS(real_datasets, save_fn)
}

# Characterize these
results <- NULL
for(i in 1:length(real_datasets)) {
  dataset_name <- datasets[i]
  data <- real_datasets[[dataset_name]]
  counts_A <- data$absolute[[1]]
  counts_B <- data$absolute[[2]]
  results <- rbind(results,
                 data.frame(dtype = dataset_name,
                            uuid = "",
                            utype = "",
                            majority_sign = get_majority_sign(counts_A, counts_B),
                            prop_overall = get_prop_sign_agreement(counts_A, counts_B),
                            p_top_10 = get_prop_top10(counts_A, counts_B),
                            max_diff = get_max_diff(counts_A, counts_B),
                            max_diff_rel = get_max_diff_relab(counts_A, counts_B),
                            p_stable = get_prop_stable(counts_A, counts_B),
                            entropy = get_entropy(counts_A, counts_B),
                            sample_agree = get_samples_agree(counts_A, counts_B),
                            zeros = get_zeros(counts_A, counts_B),
                            library_size = get_library_variation(counts_A, counts_B),
                            depth_fc = get_depth_fc(rbind(data$absolute[[1]],
                                                          data$absolute[[2]]),
                                                    rbind(data$relative[[1]],
                                                          data$relative[[2]])),
                            TP = NA,
                            TN = NA,
                            FP = NA,
                            FN = NA,
                            fpr = fprs[i]))
}

# -----------------------------------------------------------------------------
#   Simulated data
# -----------------------------------------------------------------------------

string_to_calls <- function(call_string) {
  pvals <- as.numeric(strsplit(call_string, ";")[[1]])
  pvals <- p.adjust(pvals, method = "BH")
  pvals < 0.05
}

save_fn <- "bad_good_uuids.txt"
if(file.exists(save_fn)) {
  uuid_obj <- read.table(save_fn)
  uuids <- uuid_obj$uuid
  utype <- uuid_obj$type
  uuid_string <- paste0("'", paste(uuids, collapse = "', '"), "'")
  
  conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
  res <- dbGetQuery(conn, paste0("SELECT ",
                                 "datasets.UUID AS UUID, ",
                                 "METHOD, ",
                                 "PARTIAL_INFO, ",
                                 "BASELINE_TYPE, ",
                                 "datasets.BASELINE_CALLS AS ORACLE_BASELINE, ",
                                 "CALLS, ",
                                 "results.BASELINE_CALLS AS SELF_BASELINE, ",
                                 "P, ",
                                 "CORRP, ",
                                 "LOG_MEAN, ",
                                 "PERTURBATION, ",
                                 "REP_NOISE, ",
                                 "FC_ABSOLUTE, ",
                                 "FC_RELATIVE, ",
                                 "FC_PARTIAL, ",
                                 "MED_ABS_TOTAL, ",
                                 "MED_REL_TOTAL, ",
                                 "PERCENT_DIFF_REALIZ, ",
                                 "TPR, ",
                                 "FPR ",
                                 "FROM results LEFT JOIN datasets ON ",
                                 "results.UUID=datasets.UUID ",
                                 "WHERE METHOD='DESeq2' AND PARTIAL_INFO=0 ",
                                 "AND BASELINE_TYPE='self' AND results.UUID IN (", uuid_string, ")"))
  dbDisconnect(conn)
  
  n <- nrow(res)
  tps <- numeric(n)
  tns <- numeric(n)
  fps <- numeric(n)
  fns <- numeric(n)
  for(i in 1:n) {
    baseline <- string_to_calls(res$SELF_BASELINE[i])
    calls <- string_to_calls(res$CALLS[i])
    tps[i] <- sum(baseline & calls)/length(baseline)
    tns[i] <- sum(!baseline & !calls)/length(baseline)
    fps[i] <- sum(!baseline & calls)/length(baseline)
    fns[i] <- sum(baseline & !calls)/length(baseline)
  }
  fprs <- res$FPR
} else {
  conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
  
  p <- 1000
  partial <- 0
  reference <- "self"
  
  cat(paste0("Evaluating P=", p, ", PARTIAL=", partial, ", REF=", reference, "\n"))
  res <- dbGetQuery(conn, paste0("SELECT ",
                                 "datasets.UUID AS UUID, ",
                                 "METHOD, ",
                                 "PARTIAL_INFO, ",
                                 "BASELINE_TYPE, ",
                                 "datasets.BASELINE_CALLS AS ORACLE_BASELINE, ",
                                 "CALLS, ",
                                 "results.BASELINE_CALLS AS SELF_BASELINE, ",
                                 "P, ",
                                 "CORRP, ",
                                 "LOG_MEAN, ",
                                 "PERTURBATION, ",
                                 "REP_NOISE, ",
                                 "FC_ABSOLUTE, ",
                                 "FC_RELATIVE, ",
                                 "FC_PARTIAL, ",
                                 "MED_ABS_TOTAL, ",
                                 "MED_REL_TOTAL, ",
                                 "PERCENT_DIFF_REALIZ, ",
                                 "TPR, ",
                                 "FPR ",
                                 "FROM results LEFT JOIN datasets ON ",
                                 "results.UUID=datasets.UUID ",
                                 "WHERE P=", p, " ",
                                 "AND PARTIAL_INFO=", partial, " ",
                                 "AND BASELINE_TYPE='",reference,"' ",
                                 "AND FC_ABSOLUTE <= 10 ",
                                 "AND FC_ABSOLUTE >= 0.1;"))
  
  dbDisconnect(conn)
  
  # Strip "result-less" entries
  res <- res %>%
    filter(!is.na(TPR) & !is.na(FPR))
  
  res$FC_plot <- sapply(res$FC_ABSOLUTE, function(x) {
    if(x < 1) {
      1 / x
    } else {
      x
    }
  })
  res$FC_plot <- cut(res$FC_plot, breaks = c(1, 2, 5, Inf))
  levels(res$FC_plot) <- c("low", "moderate", "high")
  
  method <- "DESeq2"; fpr_high <- 0.3; fpr_low <- 0.05
  # method <- "ALDEx2"; fpr_high <- 0.20; fpr_low <- 0.04
  # method <- "scran"; fpr_high <- 0.15; fpr_low <- 0.03
  
  # Select as "bad" data sets, those with high FC, high TPR, and high FPR
  bad_sets <- res %>%
    filter(METHOD == method) %>%
    filter(FC_plot == "high") %>%
    filter(FPR > fpr_high) %>%
    filter(TPR > 0.8) #%>%
    #select(-c("ORACLE_BASELINE", "CALLS"))
  
  # Select as "good" data sets, those with high FC, high TPR, and low FPR
  good_sets <- res %>%
    filter(METHOD == method) %>%
    filter(FC_plot == "high") %>%
    filter(FPR < fpr_low) %>%
    filter(TPR > 0.8) #%>%
    #select(-c("ORACLE_BASELINE", "CALLS"))
  
  # Do we have a reasonable sample of each?
  dim(bad_sets)
  dim(good_sets)
  
  # Is correlation any different between "bad" and "good" datasets?
  # table(bad_sets$CORRP)
  # table(good_sets$CORRP)

  uuids <- c(bad_sets$UUID, good_sets$UUID)
  fprs <- c(bad_sets$FPR, good_sets$FPR)
  utype <- c(rep("bad", nrow(bad_sets)), rep("good", nrow(good_sets)))
  
  bad_baselines <- t(sapply(bad_sets$SELF_BASELINE, function(x) string_to_calls(x)))
  rownames(bad_baselines) <- NULL
  bad_calls <- t(sapply(bad_sets$CALLS, function(x) string_to_calls(x)))
  rownames(bad_calls) <- NULL
  
  bad_tps <- numeric(nrow(bad_baselines))
  bad_tns <- numeric(nrow(bad_baselines))
  bad_fps <- numeric(nrow(bad_baselines))
  bad_fns <- numeric(nrow(bad_baselines))
  for(j in 1:nrow(bad_baselines)) {
    bad_tps[j] <- sum(bad_baselines[j,] & bad_calls[j,])
    bad_tns[j] <- sum(!bad_baselines[j,] & !bad_calls[j,])
    bad_fps[j] <- sum(!bad_baselines[j,] & bad_calls[j,])
    bad_fns[j] <- sum(bad_baselines[j,] & !bad_calls[j,])
  }
  
  good_baselines <- t(sapply(good_sets$SELF_BASELINE, function(x) string_to_calls(x)))
  rownames(good_baselines) <- NULL
  good_calls <- t(sapply(good_sets$CALLS, function(x) string_to_calls(x)))
  rownames(good_calls) <- NULL
  
  good_tps <- numeric(nrow(good_baselines))
  good_tns <- numeric(nrow(good_baselines))
  good_fps <- numeric(nrow(good_baselines))
  good_fns <- numeric(nrow(good_baselines))
  for(j in 1:nrow(good_baselines)) {
    good_tps[j] <- sum(good_baselines[j,] & good_calls[j,])
    good_tns[j] <- sum(!good_baselines[j,] & !good_calls[j,])
    good_fps[j] <- sum(!good_baselines[j,] & good_calls[j,])
    good_fns[j] <- sum(good_baselines[j,] & !good_calls[j,])
  }
  
  tps <- c(bad_tps, good_tps)
  tns <- c(bad_tns, good_tns)
  fps <- c(bad_fps, good_fps)
  fns <- c(bad_fns, good_fns)
  
  # bad_sets <- bad_sets %>%
  #   select(-c("ORACLE_BASELINE", "CALLS"))
  # good_sets <- good_sets %>%
  #   select(-c("ORACLE_BASELINE", "CALLS"))
}

output_fn <- "bad_good_results2.tsv"

# Initialize output
output_file <- file(output_fn)
writeLines(paste0(c("dtype", "uuid", "utype", "majority_sign",
                    "prop_overall", "p_top_10", "max_diff", "max_diff_rel",
                    "p_stable", "entropy", "sample_agree", "zeros",
                    "library_size", "depth_fc", "TP", "TN", "FP", "FN", "fpr"),
                  collapse = "\t"),
           output_file)
for(i in 1:nrow(results)) {
  write_delim(results[i,], output_fn, delim = "\t", append = TRUE)
}
close(output_file)

for(i in 1:length(uuids)) {
  cat(paste0("Iteration ", i, " / ", length(uuids), "\n"))
  ref_obj <- readRDS(file.path("output", "datasets", paste0(uuids[i], ".rds")))$simulation
  ref_data <- ref_obj$abundances
  groups <- ref_obj$groups
  n <- floor(nrow(ref_data)/2)
  counts_A <- ref_data[1:n,]
  counts_B <- ref_data[(n+1):(n*2),]

  results <- data.frame(dtype = "simulated",
                        uuid = uuids[i],
                        utype = utype[i],
                        majority_sign = get_majority_sign(counts_A, counts_B),
                        prop_overall = get_prop_sign_agreement(counts_A, counts_B),
                        p_top_10 = get_prop_top10(counts_A, counts_B),
                        max_diff = get_max_diff(counts_A, counts_B),
                        max_diff_rel = get_max_diff_relab(counts_A, counts_B),
                        p_stable = get_prop_stable(counts_A, counts_B),
                        entropy = get_entropy(counts_A, counts_B),
                        sample_agree = get_samples_agree(counts_A, counts_B),
                        zeros = get_zeros(counts_A, counts_B),
                        library_size = get_library_variation(counts_A, counts_B),
                        depth_fc = get_depth_fc(ref_obj$abundances,
                                                ref_obj$observed_counts1),
                        TP = tps[i],
                        TN = tns[i],
                        FP = fps[i],
                        FN = fns[i],
                        fpr = fprs[i])
  write_delim(results, output_fn, delim = "\t", append = TRUE)
}

# cat("Saving results...\n")
# saveRDS(results, "bad_good_results.rds")
