library(tidyverse)
library(codaDE)

# ------------------------------------------------------------------------------
#   Parse Athanasiadou et al. data
# ------------------------------------------------------------------------------

file_dir <- file.path("data", "Athanasiadou_2021", "S1CodeandData")
RNA_file <- "YmatRNAciona.txt"
groups <- as.factor(rep(c("lacz", "dnfgfr", "camras"), rep(3,3)))

relative_counts <- as.matrix(read.table(file = file.path(file_dir, RNA_file)))
absolute_counts <- normalize_Athanasiadou(file_dir,
                                 SI_file = "YmatSIciona.txt",
                                 RNA_file = RNA_file,
                                 ERCC_annot = "ERCCciona.txt",
                                 groups = groups)
control_group <- "lacz"
treatment_group <- "dnfgfr"
relative_counts <- relative_counts[,groups %in% c(control_group, treatment_group)]
absolute_counts <- absolute_counts[,groups %in% c(control_group, treatment_group)]

# Scale each set of counts such that minimum observed count is 1
min_observed <- min(relative_counts[relative_counts != 0])
relative_counts <- relative_counts/min_observed
min_observed <- min(absolute_counts[absolute_counts != 0])
absolute_counts <- absolute_counts/min_observed

# Round all; some methods don't tolerate floats
relative_counts <- ceiling(relative_counts)
absolute_counts <- ceiling(absolute_counts)

# Remove features that fully drop out
dropouts <- unname(which(rowSums(absolute_counts) == 0 |
                           rowSums(relative_counts) == 0))
retain_idx <- setdiff(1:nrow(relative_counts), dropouts)
relative_counts <- relative_counts[retain_idx,]
absolute_counts <- absolute_counts[retain_idx,]

cat("Approx. fold change in totals:",
    round(mean(rowSums(absolute_counts[groups == treatment_group,])) /
            mean(rowSums(absolute_counts[groups == control_group,])), 2), "\n")

# Flip to give samples x features
relative_counts <- t(relative_counts)
absolute_counts <- t(absolute_counts)

# ------------------------------------------------------------------------------
#   Apply NB GLM to absolute counts to get baseline DE
#
#   Need to think carefully about how to implement clustering in scran on real
#   data, which is much noisier than the simulations. (TBD)
# ------------------------------------------------------------------------------

# Testing.
relative_counts <- relative_counts[,1:1000]
absolute_counts <- absolute_counts[,1:1000]

groups_binary <- c(rep(0, 3), rep(1, 3))
oracle_calls <- sapply(1:ncol(absolute_counts), function(idx) {
  call_DA_NB(absolute_counts, groups_binary, idx)
})
oracle_calls <- oracle_calls < 0.05

call_DA_validation <- function(absolute_counts, relative_counts) {
  # NB GLM
  NB_calls <- sapply(1:ncol(relative_counts), function(idx) {
    call_DA_NB(relative_counts, groups_binary, idx)
  })
  # DESeq2
  DESeq2_calls <- call_DA_Seurat(relative_counts, groups_binary, method = "DESeq2")
  # MAST
  MAST_calls <- call_DA_MAST(relative_counts, groups_binary)
  # scran
  # scran_calls <- call_DA_scran(relative_counts, groups_binary, do_cluster = FALSE)
  # ALDEx2
  ALDEx2_calls <- call_DA_ALDEx2(relative_counts, groups_binary)
  calls_list <- list(baseline = NB_calls,
                     DESeq2 = DESeq2_calls,
                     MAST = MAST_calls,
                     # scran = scran_calls,
                     ALDEx2 = ALDEx2_calls)
  return(calls_list)
}

calls_list <- call_DA_validation(absolute_counts, relative_counts)

calculate_rates_validation <- function(calls_list) {
  results <- data.frame(tpr = c(),
                        fpr = c(),
                        method = c())
  for(method_idx in 1:length(calls_list)) {
    DE_calls <- calls_list[[method_idx]]
    DE_calls <- DE_calls < 0.05
    
    TP <- sum(oracle_calls & DE_calls)
    FP <- sum(!oracle_calls & DE_calls)
    
    TN <- sum(!oracle_calls & !DE_calls)
    FN <- sum(oracle_calls & !DE_calls)
    
    TPR <- TP / (TP + FN)
    FPR <- FP / (FP + TN)
    
    results <- rbind(results,
                     data.frame(tpr = TPR,
                                fpr = FPR,
                                method = names(calls_list)[method_idx]))
  }
  return(results)
}

results <- calculate_rates_validation(calls_list)
saveRDS(results, file = file.path("output", "Athanasiadou_validation_results.rds"))

