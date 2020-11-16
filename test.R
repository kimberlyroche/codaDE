library(codaDE)

p <- 1000
n <- 100
prop_da <- 0.9
k <- 1
sf_corr <- 0
sim <- simulate_singlecell_RNAseq(p = p, n = n, k = 1, proportion_da = prop_da, size_factor_correlation = sf_corr, spike_in = FALSE, possible_fold_changes = c(6, 8, 10, 50))

# Using the NB with no correction give huge error
data <- evaluate_DA(sim, alpha = 0.05, use_ALR = FALSE, filter_abundance = 0, method = "NB")
data$FPR
data$TPR
plot(data$library_size.abundances)

# Using edgeR with dispersion estimation, etc. hugely reduces error
data <- evaluate_DA(sim, alpha = 0.05, use_ALR = FALSE, filter_abundance = 0, method = "edgeR")
data$FPR
data$TPR
plot(data$library_size.abundances)
