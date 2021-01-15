library(ggplot2)
library(gridExtra)

GTEx_estimates_pairs <- readRDS("data/empirical_mu_size_pairs.rds")

idx <- sample(1:nrow(GTEx_estimates_pairs), size = 1)
y1 <- rnbinom(100, mu = GTEx_estimates_pairs[idx,1], size = GTEx_estimates_pairs[idx,2])
y2 <- rnbinom(100, mu = GTEx_estimates_pairs[idx,3], size = GTEx_estimates_pairs[idx,4])

ggplot(data.frame(y = c(y1, y2), x = 1:200, tissue = as.factor(c(rep("A", 100), rep("B", 100)))), aes(x = x, y = y, color = tissue)) +
  geom_point(size = 2) +
  xlab("sample index") +
  ylab("abundance")

library(codaDE)

# Simulate from a fixed fold change library
# Note: sequencing depth with determine the number of genes you have *resolution* to detect
# Number of simulated genes is 20K; number of non-zero genes is typically much less
# sim_data <- simulate_singlecell_RNAseq(n = 200, sequencing_depth = 1e4, proportion_da = 0.7, library_size_correlation = 0,
#                                        spike_in = FALSE, possible_fold_changes = c(2, 4, 8))

# Simulate from empirical data (w/ cutoff)
sim_data <- simulate_singlecell_RNAseq(n = 200, sequencing_depth = 1e4, proportion_da = 0.7, library_size_correlation = 0,
                                   spike_in = FALSE)

# Takes ~2min. with sequencing depth = 1e5
calls <- call_DA_edgeR(sim_data, call_abundances = FALSE, normalization_method = NULL)
str(calls)

de <- which(calls < 0.05)
not_de <- setdiff(1:length(calls), de)

FP <- setdiff(de, sim_data$da_assignment)
FN <- setdiff(not_de, sim_data$da_assignment)

# FPR
length(FP) / 20000

# FNR
length(FN) / 20000

# Delta total abundance
plot(rowSums(sim_data$abundances))
plot(rowSums(sim_data$observed_counts))

idx <- sample(FP, size = 1)
p1 <- ggplot(data.frame(y = c(sim_data$abundances[,idx]),
                        x = 1:400,
                        tissue = as.factor(c(rep("A", 200), rep("B", 200)))), aes(x = x, y = y, color = tissue)) +
  geom_point(size = 2) +
  xlab("sample index") +
  ylab("abundance")
p2 <- ggplot(data.frame(y = c(sim_data$observed_counts[,idx]),
                        x = 1:400,
                        tissue = as.factor(c(rep("A", 200), rep("B", 200)))), aes(x = x, y = y, color = tissue)) +
  geom_point(size = 2) +
  xlab("sample index") +
  ylab("abundance")
grid.arrange(grobs = list(p1, p2), ncol = 2)
