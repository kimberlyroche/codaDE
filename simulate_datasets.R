# This scripts simulates data sets under the regime with the most variation in error.
#
# Data sets with (1) lots of differentially abundant components and (2) no relationship between
# the observed and true absolute abundances can have low-ish error or lots of error.
#
# I'm simulating some data sets to see if I can discern any obvious characteristics, produced
# under this regime, that distinguish the OK from the awful.
#
# 8/18/2020

library(codaDE)

p <- 200

for(j in 1:20) {

  data <- simulate_singlecell_RNAseq(p = p, n = 250, k = 1, proportion_da = 0.8,
                                     size_factor_correlation = 0, spike_in = FALSE)

  calls.abundances <- c()
  calls.observed_counts <- c()
  alpha <- 0.05
  for(i in 1:p) {
    if(i %% 100 == 0) {
      cat("Evaluating DA on feature:",i,"\n")
    }
    pval.abundances <- call_DA(data, i, call_abundances = TRUE)
    pval.observed_counts <- call_DA(data, i, call_abundances = FALSE)
    if(!is.na(pval.abundances) & !is.na(pval.observed_counts)) {
      calls.abundances <- c(calls.abundances, pval.abundances <= alpha/p)
      calls.observed_counts <- c(calls.observed_counts, pval.observed_counts <= alpha/p)
    }
  }

  FN <- sum(calls.abundances == TRUE & calls.observed_counts == FALSE)
  FP <- sum(calls.abundances == FALSE & calls.observed_counts == TRUE)
  TN <- sum(calls.abundances == FALSE & calls.observed_counts == FALSE)
  TP <- sum(calls.abundances == TRUE & calls.observed_counts == TRUE)

  saveRDS(list(data = data, FN = FN, FP = FP, TN = TN, TP = TP,
          calls.abundances = calls.abundances, calls.observed_counts = calls.observed_counts),
          file = paste0("sim_dataset_",j,".rds"))

}
