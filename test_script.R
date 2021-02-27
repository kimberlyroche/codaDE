# 2021-02-27
# This is a constantly-edited test script for generating simulated data sets and evaluating their properties.

library(tidyverse)
library(gridExtra)
library(codaDE)

# ------------------------------------------------------------------------------------------------------------
#   Simulate data
# ------------------------------------------------------------------------------------------------------------

n <- 25
p <- 200
ref_data <- "test"
k <- 1
sequencing_depth <- 1e5
proportion_da <- 0.5
library_size_correlation <- 0
spike_in <- FALSE
possible_fold_changes <- NULL

sim_data <- simulate_sequence_counts(n = n, p = p, k = k, ref_data = ref_data,
                                     sequencing_depth = sequencing_depth, proportion_da = proportion_da,
                                     library_size_correlation = library_size_correlation,
                                     spike_in = spike_in, possible_fold_changes = possible_fold_changes)

# ------------------------------------------------------------------------------------------------------------
#   Visualize simulation - relative abundances
# ------------------------------------------------------------------------------------------------------------

# `sim_data$abudances` is (2*n) x p (samples x features)
props <- t(apply(sim_data$abundances, 1, function(x) x/sum(x)))

palette <- generate_highcontrast_palette(ncol(props))

data <- pivot_longer(cbind(sample = 1:nrow(props), as.data.frame(props)),
                     cols = !sample,
                     names_to = "OTU",
                     values_to = "relative_abundance")
ggplot(data, aes(fill = OTU, y = relative_abundance, x = sample)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = palette) +
  theme(legend.position = "none")

# ------------------------------------------------------------------------------------------------------------
#   Visualize simulation - absolute abundances
# ------------------------------------------------------------------------------------------------------------

data <- pivot_longer(cbind(sample = 1:nrow(sim_data$abundances), as.data.frame(sim_data$abundances)),
                     cols = !sample,
                     names_to = "OTU",
                     values_to = "abundance")

palette <- generate_highcontrast_palette(length(unique(data$OTU)))

ggplot(data, aes(fill = OTU, y = abundance, x = sample)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = palette) +
  theme(legend.position = "none")

# ------------------------------------------------------------------------------------------------------------
#   Call differential expression/abundance
# ------------------------------------------------------------------------------------------------------------

if(FALSE) {
  calls <- call_DA_edgeR(sim_data, call_abundances = FALSE, normalization_method = NULL)

  de <- which(calls < 0.05)
  not_de <- setdiff(1:length(calls), de)

  FP <- setdiff(de, sim_data$da_assignment)
  FN <- setdiff(not_de, sim_data$da_assignment)

  # FPR
  cat("FP proportion:", round(length(FP) / ncol(sim_data$abundances), 2), "\n")

  # FNR
  cat("FN proportion:", round(length(FN) / ncol(sim_data$abundances), 2), "\n")
}
