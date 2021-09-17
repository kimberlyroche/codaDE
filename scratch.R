source("path_fix.R")

library(codaDE)
library(tidyverse)

percent_stable <- function(data, groups, abs_fc_threshold = 1.5) {
  means_A <- colMeans(data[groups == levels(groups)[1],])
  means_B <- colMeans(data[groups == levels(groups)[2],])
  fc <- (means_B + 0.1) / (means_A + 0.1)
  fc <- sapply(fc, function(x) ifelse(x < 1, 1/x, x))
  sum(fc < abs_fc_threshold)/length(fc)
}

statistics <- NULL
datasets <- c("VieiraSilva", "Barlow") #, "Song", "Monaco",
	          # "Hagai", "Owens", "Klein", "Yu")
for(dataset in datasets) {
  for(threshold in 0:5) {
  	abs_data <- do.call(paste0("parse_", dataset_name), list(absolute = TRUE))
  	rel_data <- do.call(paste0("parse_", dataset_name), list(absolute = FALSE))
  
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
  
  	n_features <- ncol(ref_data)
  	zeros_abs <- sum(ref_data == 0)/(nrow(ref_data)*ncol(ref_data))
  	zeros_rel <- sum(data == 0)/(nrow(data)*ncol(data))
  	stable_abs <- percent_stable(ref_data, groups, abs_fc_threshold = 1.5)
  	stable_rel <- percent_stable(data, groups, abs_fc_threshold = 1.5)
  	
  	statistics <- rbind(statistics,
  	                    data.frame(dataset = dataset,
  	                               filter = threshold,
  	                               n_features = n_features,
  	                               zeros_abs = zeros_abs,
  	                               zeros_rel = zeros_rel,
  	                               stable_abs = stable_abs,
  	                               stable_rel = stable_rel))
  }
}
saveRDS(statistics, file.path("output", "real_data_filtering_statistics.rds"))

