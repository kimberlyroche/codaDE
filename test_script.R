# 2021-02-27
# This is a constantly-edited test script for generating simulated data sets and evaluating their properties.

library(tidyverse)
library(gridExtra)
library(codaDE)

plot_stacked_bars <- function(data, palette, save_name = NULL) {
  ggplot(data, aes(fill = OTU, y = abundance, x = sample)) + 
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = palette) +
    theme(legend.position = "none")
  if(!is.null(save_name)) {
    ggsave(file.path("images", save_name),
           units = "in",
           dpi = 100,
           height = 5,
           width = 6)
  }
}

# ------------------------------------------------------------------------------------------------------------
#   Simulate data
# ------------------------------------------------------------------------------------------------------------

build_simulated_reference(p = 1000, log_var = 2, log_noise_var = 2)

n <- 5 # replicate number
p <- 100
palette <- generate_highcontrast_palette(p)
# ref_data <- "simulated"
# ref_data <- "Morton"
ref_data <- "Athanasiadou_ciona"
k <- 1
asymmetry <- 0.5
sequencing_depth <- 1e5
proportion_da <- 0.5
library_size_correlation <- 0
spike_in <- FALSE
possible_fold_changes <- NULL

q <- 1000
fcs <- numeric(q)
for(i in 1:q) {
  sim_data <- simulate_sequence_counts(n = n, p = p, k = k, ref_data = ref_data, asymmetry = asymmetry,
                                       sequencing_depth = sequencing_depth, proportion_da = proportion_da,
                                       library_size_correlation = library_size_correlation,
                                       spike_in = spike_in, possible_fold_changes = possible_fold_changes)
  totals <- rowSums(sim_data$abundances)
  t1 <- mean(totals[1:n])
  t2 <- mean(totals[(n+1):(n*2)])
  fc[i] <- t2 / t1
}
hist(q)

# ------------------------------------------------------------------------------------------------------------
#   Visualize simulation - absolute abundances
# ------------------------------------------------------------------------------------------------------------

data <- pivot_longer(cbind(sample = 1:nrow(sim_data$abundances), as.data.frame(sim_data$abundances)),
                     cols = !sample,
                     names_to = "OTU",
                     values_to = "abundance")

plot_stacked_bars(data, palette, save_name = paste0("01_absolute_",ref_data,".png"))

# ------------------------------------------------------------------------------------------------------------
#   Visualize simulation - observed abundances
# ------------------------------------------------------------------------------------------------------------

data <- pivot_longer(cbind(sample = 1:nrow(sim_data$observed_counts), as.data.frame(sim_data$observed_counts)),
                     cols = !sample,
                     names_to = "OTU",
                     values_to = "abundance")

plot_stacked_bars(data, palette, save_name = paste0("02_observed_",ref_data,".png"))

# ------------------------------------------------------------------------------------------------------------
#   Visualize simulation - relative abundances
# ------------------------------------------------------------------------------------------------------------

# `sim_data$abundances` is (2*n) x p (samples x features)
props <- t(apply(sim_data$abundances, 1, function(x) x/sum(x)))
data <- pivot_longer(cbind(sample = 1:nrow(props), as.data.frame(props)),
                     cols = !sample,
                     names_to = "OTU",
                     values_to = "abundance") # really, relative abundance

plot_stacked_bars(data, palette, save_name = paste0("03_relative_",ref_data,".png"))

# ------------------------------------------------------------------------------------------------------------
#   Call differential expression/abundance
# ------------------------------------------------------------------------------------------------------------

oracle_calls <- sapply(1:p, function(idx) { call_DA_NB(sim_data, idx, call_abundances = TRUE) } )
calls <- sapply(1:p, function(idx) { call_DA_NB(sim_data, idx, call_abundances = FALSE) } )

de <- calls < 0.05
sim_de <- oracle_calls < 0.05

TP <- sum(de & sim_de)
FP <- sum(de & !sim_de)

TN <- sum(!de & !sim_de)
FN <- sum(!de & sim_de)

# FPR
cat("TPR:", round(TP/(TP+FN), 2), "\n")

# FNR
cat("FPR:", round(FP/(FP+TN), 2), "\n")
