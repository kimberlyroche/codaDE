source("path_fix.R")

library(codaDE)

datasets_names <- c("VieiraSilva", "Barlow", "Song", "Monaco", "Hagai", "Owens", "Klein", "Yu")
thresholds <- c(1, 1, 1, 2, 1, 1, 1, 1)

for(i in 1:length(datasets_names)) {
  name <- datasets_names[i]
  ts <- thresholds[i]

  abs_data <- do.call(paste0("parse_", name), list(absolute = TRUE))
  rel_data <- do.call(paste0("parse_", name), list(absolute = FALSE))
  
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
  retain_features <- colMeans(ref_data) >= ts & colMeans(data) >= ts
  ref_data <- ref_data[,retain_features]
  data <- data[,retain_features]
  
  cat(paste0("Dataset dimensions: ", nrow(ref_data), " samples x ", ncol(ref_data), " features\n"))

  fc <- mean(rowSums(ref_data[groups == unique(groups)[1],])) / mean(rowSums(ref_data[groups == unique(groups)[2],]))
  if(fc < 1) {
    fc <- 1 / fc
  }

  cat(paste0("\tFC: ", round(fc, 1), "\n"))
  cat(paste0("\tPercent DE (thresholded): ", 100*round(sum(calc_threshold_DA(ref_data, nA = sum(groups == unique(groups)[1]))) / ncol(ref_data), 2), "\n"))

}
