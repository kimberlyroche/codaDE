library(codaDE)
library(tidyverse)

p <- 100
n <- 50

iter <- 100
asymmetry <- 1
proportion_da <- 0.95
spike_in <- FALSE

A <- matrix(0.8, 100, 100)
diag(A) <- 1
B <- diag(100) # indep

change_A <- c()
change_B <- c()
for(i in 1:iter) {
  # sim data
  build_simulated_reference(p = p, log_mean = 0, log_var = 2, log_noise_var = 1,
                            base_correlation = A, concentration = p + 100)
  # ref_file <- file.path("data", paste0("DE_reference_simulated.rds"))
  # data_obj <- readRDS(ref_file)

  sim_data <- simulate_sequence_counts(n = n, p = p, ref_data = "simulated",
                                       asymmetry = asymmetry,
                                       proportion_da = proportion_da,
                                       spike_in = spike_in)
  
  # change_A <- c(change_A, sum(data_obj$log_perturbation))
  # change_A <- c(change_A, sum(mean(data_obj$log_cond1) - mean(data_obj$log_cond2)))
  # change_A <- c(change_A, sum(data_obj$cond2))
  # change_A <- c(change_A, sum(sim_data$abundances[(n+1):(n*2),1:p]))
  # change_A <- c(change_A, abs(abs(sum(sim_data$baseline_counts) - sum(sim_data$treatment_counts))))
  change_A <- c(change_A, abs(sum(sim_data$abundances[1:n,1:p]) - sum(sim_data$abundances[(n+1):(n*2),1:p])))
  
  build_simulated_reference(p = p, log_mean = 0, log_var = 2, log_noise_var = 1,
                            base_correlation = B, concentration = p + 100)
  # ref_file <- file.path("data", paste0("DE_reference_simulated.rds"))
  # data_obj <- readRDS(ref_file)
  
  sim_data <- simulate_sequence_counts(n = n, p = p, ref_data = "simulated",
                                       asymmetry = asymmetry,
                                       proportion_da = proportion_da,
                                       spike_in = spike_in)
  
  # change_B <- c(change_B, sum(data_obj$log_perturbation))
  # change_B <- c(change_B, sum(mean(data_obj$log_cond1) - mean(data_obj$log_cond2)))
  # change_B <- c(change_B, sum(data_obj$cond2))
  # change_B <- c(change_B, sum(sim_data$abundances[(n+1):(n*2),1:p]))
  # change_B <- c(change_B, abs(abs(sum(sim_data$baseline_counts) - sum(sim_data$treatment_counts))))
  change_B <- c(change_B, abs(sum(sim_data$abundances[1:n,1:p]) - sum(sim_data$abundances[(n+1):(n*2),1:p])))
}

ggplot(data.frame(x = c(change_A, change_B),
                  type = c(rep("A", iter), rep("B", iter))), aes(x = x, fill = type)) +
  geom_density(alpha = 0.5) #+
  #xlim(c(0, 20000))

mean(change_A)
mean(change_B)

sd(change_A)
sd(change_B)

# In general the effect of correlation on the variance in total abundances is
# pretty small.





