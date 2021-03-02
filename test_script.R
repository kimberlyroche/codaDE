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

calc_DE_discrepancy <- function(ref_data, data, groups) {
  oracle_calls <- sapply(1:p, function(idx) { call_DA_NB(ref_data, groups, idx) } )
  calls <- sapply(1:p, function(idx) { call_DA_NB(data, groups, idx) } )

  de <- calls < 0.05
  sim_de <- oracle_calls < 0.05
  
  TP <- sum(de & sim_de)
  FP <- sum(de & !sim_de)
  
  TN <- sum(!de & !sim_de)
  FN <- sum(!de & sim_de)

  return(list(tpr = TP/(TP+FN), fpr = FP/(FP+TN)))  
}

# ------------------------------------------------------------------------------------------------------------
#   Simulate data
# ------------------------------------------------------------------------------------------------------------

# build_simulated_reference(p = 1000, log_var = 2, log_noise_var = 2)

n <- 5 # replicate number
p <- 100
palette <- generate_highcontrast_palette(p)
ref_data <- "simulated"
# ref_data <- "Morton"
# ref_data <- "Barlow"
# ref_data <- "Athanasiadou_ciona"
# ref_data <- "Athanasiadou_yeast"
k <- 1
asymmetry <- 0.5
sequencing_depth <- 1e5
proportion_da <- 0.5
spike_in <- FALSE
possible_fold_changes <- NULL

if(FALSE) { # visualize distribution of fold changes from 1K simulations
  q <- 1000
  fcs <- numeric(q)
  for(i in 1:q) {
    if(i %% 100 == 0) {
      cat("Iteration",i,"\n")
    }
    sim_data <- simulate_sequence_counts(n = n, p = p, k = k, ref_data = ref_data, asymmetry = asymmetry,
                                         sequencing_depth = sequencing_depth, proportion_da = proportion_da,
                                         spike_in = spike_in, possible_fold_changes = possible_fold_changes)
    totals <- rowSums(sim_data$abundances)
    t1 <- mean(totals[1:n])
    t2 <- mean(totals[(n+1):(n*2)])
    fc[i] <- t2 / t1
  }
  
  data <- data.frame(fc = fc)
  ggplot(data, aes(x = fc)) +
    geom_histogram(color = "white") +
    xlim(0, 5)
  ggsave(file.path("images", paste0("00_hist_",ref_data,"_p",p,".png")),
         units = "in",
         dpi = 100,
         height = 5,
         width = 6)
} else {
  sim_data <- simulate_sequence_counts(n = n, p = p, k = k, ref_data = ref_data, asymmetry = asymmetry,
                                       sequencing_depth = sequencing_depth, proportion_da = proportion_da,
                                       spike_in = spike_in, possible_fold_changes = possible_fold_changes)
}

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

# Scenario 1: No correlation between true abundances and observed abundances
data <- pivot_longer(cbind(sample = 1:nrow(sim_data$observed_counts1), as.data.frame(sim_data$observed_counts1)),
                     cols = !sample,
                     names_to = "OTU",
                     values_to = "abundance")

plot_stacked_bars(data, palette, save_name = paste0("02_observed_",ref_data,".png"))

# Scenario 2: Partial correlation between true abundances and observed abundances
data <- pivot_longer(cbind(sample = 1:nrow(sim_data$observed_counts2), as.data.frame(sim_data$observed_counts2)),
                     cols = !sample,
                     names_to = "OTU",
                     values_to = "abundance")

plot_stacked_bars(data, palette, save_name = paste0("02_observed_",ref_data,"_partialcorr.png"))

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

rates1 <- calc_DE_discrepancy(sim_data$abundances, sim_data$observed_counts1, sim_data$groups)
rates2 <- calc_DE_discrepancy(sim_data$abundances, sim_data$observed_counts2, sim_data$groups)

rates1
rates2
