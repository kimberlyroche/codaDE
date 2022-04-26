source("path_fix.R")

library(codaDE)
library(ggplot2)
library(cowplot)

threshold <- 1
model_subdir <- "regression_plain"
use_cpm <- FALSE

datasets <- c("Song", "Monaco", "Muraro", "Hagai", "Hashimshony", "Gruen")
for(dataset_name in datasets) {
  parsed_obj <- wrangle_validation_data(dataset_name = dataset_name,
                                        threshold = threshold,
                                        use_cpm = use_cpm,
                                        testing = FALSE,
                                        hkg_list = hkg)
  ref_data <- parsed_obj$ref_data
  data <- parsed_obj$data
  groups <- parsed_obj$groups
  tax <- parsed_obj$tax
  hkg_present <- parsed_obj$hkg_present
  hkg_rho <- parsed_obj$hkg_rho
  
  calls <- readRDS(paste0("output/real_data_calls/no_norm/calls_oracle_DESeq2_",dataset_name,"_threshold1.rds"))
  calls_cg <- readRDS(paste0("output/real_data_calls/no_norm/calls_oracle_DESeq2_CG_",dataset_name,"_threshold1_noHKG.rds"))
  
  cat(paste0("Dataset: ", dataset_name, "\n"))
  cat(paste0("\tRho: ", round(hkg_rho, 3), "\n"))
  cat(paste0("\tTPR: ",
             round(calls$rates$TPR, 3), " -> ",
             round(calls_cg$rates$TPR, 3), " (FC = ",
             round(calls_cg$rates$TPR/calls$rates$TPR, 3), ")\n"))
  cat(paste0("\tFPR: ", round(calls$rates$FPR, 3), " -> ", round(calls_cg$rates$FPR, 3), " (FC = ",
             round(calls_cg$rates$FPR/calls$rates$FPR, 3), ")\n"))
}

# Get these from `call_DA_DESeq2`
# dds_data_before <- t(counts(dds, normalized = TRUE))
# sf_before <- sizeFactors(dds)
# dds_data_after <- t(counts(dds, normalized = TRUE))
# sf_after <- sizeFactors(dds)

# What are the Song et al. misses that are recovered?
# These are borderline
fixed <- which(calls$rates$FN_calls & calls_cg$rates$TP_calls) # better
fixed <- which(calls$rates$TN_calls & calls_cg$rates$FP_calls) # worse
idx <- sample(fixed, size = 1)
p1 <- ggplot(data.frame(x = 1:nrow(data),
                  y = dds_data_before[,idx],
                  group = groups),
       aes(x = x, y = y, fill = group)) +
  geom_point(size = 3, shape = 21) +
  labs(title = "no CG")
p2 <- ggplot(data.frame(x = 1:nrow(data),
                        y = dds_data_after[,idx],
                        group = groups),
             aes(x = x, y = y, fill = group)) +
  geom_point(size = 3, shape = 21) +
  labs(title = "with CG")
p3 <- ggplot(data.frame(x = 1:nrow(data),
                        y = ref_data[,idx],
                        group = groups),
             aes(x = x, y = y, fill = group)) +
  geom_point(size = 3, shape = 21) +
  labs(title = "truth")
plot_grid(p1, p2, p3, ncol = 3)

plot(rowSums(ref_data), 1/sf_before)
plot(rowSums(ref_data), 1/sf_after)
cor(rowSums(ref_data), 1/sf_before)
cor(rowSums(ref_data), 1/sf_after)




