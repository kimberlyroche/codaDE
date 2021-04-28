library(codaDE)
library(tidyverse)

p <- 100
n <- 10

iter <- 500
asymmetry <- 0.8
proportion_da <- 0.75
spike_in <- TRUE

mean_corr <- numeric(iter)
lfc <- numeric(iter)
for(i in 1:iter) {
  # sim data
  build_simulated_reference(p = p, log_mean = 0, log_var = 2, log_noise_var = 1,
                            base_correlation = 0.75, concentration = p + 100)
  # ref_file <- file.path("data", paste0("DE_reference_simulated.rds"))
  # data_obj <- readRDS(ref_file)
  
  sim_data <- simulate_sequence_counts(n = n, p = p, ref_data = "simulated_bulk",
                                       asymmetry = asymmetry,
                                       proportion_da = proportion_da,
                                       spike_in = spike_in)
  
  # eyeball what observed log-counts would look like
  # ggplot(data.frame(x = log(round(data_obj$cond1) + 1)), aes(x = x)) +
  #   geom_histogram()
  
  # expr1 <- log(data_obj$cond1 + 1)
  # expr2 <- log(data_obj$cond2 + 1)
  # subset_idx <- sample(1:p, size = 100)
  # plot_bipartite_graph(expr1[subset_idx], expr2[subset_idx], alpha = 0.5)

  # K <- data_obj$correlation_matrix
  # off_diag_corr <- K[upper.tri(K, diag = FALSE)]
  # mean_corr[i] <- mean(off_diag_corr)
  
  m1 <- mean(rowSums(sim_data$abundances[1:n,]))
  m2 <- mean(rowSums(sim_data$abundances[(n+1):(n*2),]))
  log_delta_mean <- log(max(c(m1, m2)) / min(c(m1, m2)))
  lfc[i] <- log_delta_mean
}

# par(mfrow = c(1,2))
# plot(mean_corr, lfc, xlim = c(0, 2))
hist(lfc, xlim = c(0, 2))

