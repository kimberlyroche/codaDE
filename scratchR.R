library(codaDE)
library(gridExtra)

p <- 200
n <- 10

half_p <- round(p/2)
base_correlation <- matrix(0, p, p)

# ------------------------------------------------------------------------------
#   Sweep over 4 combos of correlated features
# ------------------------------------------------------------------------------

# Allow correlation to have levels
rhos <- c(0.15, 0.3, 0.45, 0.6)
concentrations <- seq(from = 50, to = 20, length.out = 4)

pl <- list()
for(corrp in 1:4) {
  if(is.null(corrp)) {
    # Independent features
    base_correlation <- diag(p)
    concentration <- 1e6
  } else {
    base_correlation[1:half_p,1:half_p] <- rhos[corrp]
    diag(base_correlation) <- 1
    concentration <- p + concentrations[corrp]
  }
  
  K <- cov2cor(rInvWishart(1, concentration, base_correlation*concentration)[,,1])
  print(median(abs(K[upper.tri(K, diag = FALSE)])))
  pl[[corrp]] <- plot_kernel_or_cov_matrix(K)
  # plot(density(K[upper.tri(K, diag = FALSE)]))
}
grid.arrange(grobs = pl, ncol = 2)
plot(density(K[upper.tri(K, diag = FALSE)]))

# ------------------------------------------------------------------------------
#   Sweep over perturbations sizes
# ------------------------------------------------------------------------------

perturbations <- seq(from = 0.1, to = 4, length.out = 10)
i <- 5

# ------------------------------------------------------------------------------
#   Sweep over log mean and replicate noise
# ------------------------------------------------------------------------------

# log_mean should range from 2 to 6
# replicate_noise should range from 0 to 1

# Test extremes; make sure they're not batshit
log_mean <- 2
replicate_noise <- 1

# Create data set
data_obj <- build_simulated_reference(p = p,
                                      log_mean = log_mean,
                                      log_noise_var = perturbations[i],
                                      base_correlation = base_correlation,
                                      concentration = concentration)
sim_data <- simulate_sequence_counts(n = n,
                                     p = p,
                                     data_obj = data_obj,
                                     replicate_noise = replicate_noise)

# plot_stacked_bars(sim_data$abundances)
plot_stacked_bars(sim_data$observed_counts1)

# Percent zeros
sum(sim_data$abundances == 0)/(n*p)
sum(sim_data$observed_counts1 == 0)/(n*p)









