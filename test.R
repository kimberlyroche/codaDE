library(codaDE)

p <- 1000
n <- 100
prop_da <- 0.8
k <- 1
sf_corr <- 0
sim_dir <- "simulated_analyses"

data <- simulate_singlecell_RNAseq(p = p, n = n, k = k, proportion_da = prop_da, size_factor_correlation = sf_corr)
res <- evaluate_DA(data, alpha = 0.05, use_ALR = FALSE, filter_abundance = 0, method = "edgeR")
res$FPR
res$TPR

res <- evaluate_DA(data, alpha = 0.05, use_ALR = FALSE, filter_abundance = 0, method = "edgeR_scran")
res$FPR
res$TPR

res <- evaluate_DA(data, alpha = 0.05, use_ALR = FALSE, filter_abundance = 0, method = "edgeR_TMM")
res$FPR
res$TPR



# Sweep
sweep_simulations(p, n, k = k, de_sweep = c(prop_da), corr_sweep = c(sf_corr),
                  alpha = 0.05, use_ALR = FALSE, filter_abundance = 0, methods = c("edgeR"))

output_files <- list.files(sim_dir)
for(i in 1:2) {
  data <- readRDS(file.path(sim_dir, output_files[i]))
  cat(paste0("Method: ",data$method, ", FPR: ",round(data$FPR, 3),"\n"))
  if(i == 1) {
    plot(data$library_size.abundances)
  }
}
