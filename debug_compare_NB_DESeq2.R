library(codaDE)
library(ggplot)

calc_fc <- function(M) {
  counts_A <- M[1:(nrow(M)/2),]
  counts_B <- M[(nrow(M)/2+1):nrow(M),]
  m1 <- mean(rowSums(counts_A))
  m2 <- mean(rowSums(counts_B))
  m2 / m1
  # max(c(m1, m2)) / min(c(m1, m2))
}

P <- 1000
CORRP <- 2

half_p <- round(P/2) # need this repeatedly for later calculations

base_correlation <- matrix(0, P, P)
# Partial correlation
base_correlation[1:half_p,1:half_p] <- 0.3
diag(base_correlation) <- 1
concentration <- P + 40

fc <- 0
while(fc < 2.5) {
  cat("Simulating dataset...\n")
  
  # Create reference distributions
  data_obj <- build_simulated_reference(p = P,
                                        log_mean = 4,
                                        perturb_size = 1,
                                        base_correlation = base_correlation,
                                        concentration = concentration)
  
  # Sample from these to create the data set
  sim_data <- simulate_sequence_counts(n = 10,
                                       p = P,
                                       data_obj = data_obj,
                                       replicate_noise = 0.6,
                                       proportion_da = 0.5)
  
  fc <- calc_fc(sim_data$abundances)
}

NB_calls <- p.adjust(call_DA_NB(sim_data$abundances, sim_data$groups)$pval, method = "BH")
DESeq2_calls <- p.adjust(call_DA_DESeq2(sim_data$abundances, sim_data$groups)$pval, method = "BH")

# What does DESeq2 call that the NB GLM doesn't?
idx <- which(NB_calls >= 0.05 & DESeq2 < 0.05)
ggplot(data.frame(x = 1:20, y = sim_data$abundances[,sample(idx, size = 1)], type = factor(sim_data$groups)),
       aes(x = x, y = y, fill = type)) +
  geom_point(size = 2, shape = 21)


