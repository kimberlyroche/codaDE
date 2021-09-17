source("path_fix.R")

library(codaDE)

datasets_names <- c("VieiraSilva", "Barlow", "Song", "Monaco", "Hagai", "Owens", "Klein", "Yu")
thresholds <- c(1, 1, 1, 2, 1, 1, 1, 1)
for(i in 1:length(datasets_names)) {
  name <- datasets_names[i]
  ts <- thresholds[i]
  data <- readRDS(file.path("output", paste0("calls_self_DESeq2_", name, "_threshold", ts, ".rds")))
  n <- length(data$all_calls$oracle_calls)
  baseline_yes <- sum(p.adjust(data$all_calls$oracle_calls, method = "BH") < 0.05)/n
  tp <- sum(data$rates$TP_calls)
  fp <- sum(data$rates$FP_calls)
  tn <- sum(data$rates$TN_calls)
  fn <- sum(data$rates$FN_calls)
  cat(paste0(name, "\n\tDE %: ", round(baseline_yes, 2), "\n\tFPR: ", round(fp/(fp+tn), 3), "\n\tFP %: ", round(fp/n, 3), "\n\tTN %: ", round(tn/n, 3), "\n"))
  cat(paste0("\tFP: ", fp, "\n\tTN: ", tn, "\n"))

  # Check what NB calls look like; code copied from validate
  abs_data <- do.call(paste0("parse_", name), list(absolute = TRUE))
  rel_data <- do.call(paste0("parse_", name), list(absolute = FALSE))

  if(nrow(abs_data$counts) > 500) {
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

  A_idx <- groups == levels(groups)[1]
  B_idx <- groups == levels(groups)[2]

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
  retain_features <- colMeans(ref_data) >= ts & colMeans(data) >= ts
  ref_data <- ref_data[,retain_features]
  # data <- data[,retain_features]

  cat(paste0("Dataset dimensions: ", nrow(ref_data), " x ", ncol(ref_data), "\n"))
  # cat(paste0("Percent zeros: ", round(sum(data == 0)/(nrow(data)*ncol(data)), 3)*100, "%\n"))

  fc <- mean(rowSums(ref_data[B_idx,])) / mean(rowSums(ref_data[A_idx,]))
  if(fc < 1) {
    fc <- 1 / fc
  }
  cat(paste0("\tFC: ", round(fc, 1), "\n"))
  cat(paste0("\tPercent DE (oracle): ", 100*round(sum(p.adjust(call_DA_NB(ref_data, groups)$pval, method = "BH") < 0.05)/ncol(ref_data), 3), "\n"))
}
