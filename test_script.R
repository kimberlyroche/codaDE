# 2021-02-27
# This is a constantly-edited test script for generating simulated data sets and evaluating their properties.

library(tidyverse)
library(gridExtra)
library(codaDE)

# ------------------------------------------------------------------------------------------------------------
#   Functions
# ------------------------------------------------------------------------------------------------------------

plot_stacked_bars <- function(data, palette, save_name = NULL) {
  p <- ggplot(data, aes(fill = OTU, y = abundance, x = sample)) + 
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = palette) +
    theme(legend.position = "none")
  show(p)
  if(!is.null(save_name)) {
    ggsave(file.path("images", save_name),
           p,
           units = "in",
           dpi = 100,
           height = 5,
           width = 6)
  }
}

DE_by_NB <- function(ref_data, data, groups) {
  oracle_calls <- sapply(1:p, function(idx) { call_DA_NB(ref_data, groups, idx) } )
  calls <- sapply(1:p, function(idx) { call_DA_NB(data, groups, idx) } )
  return(list(oracle_calls = oracle_calls, calls = calls))
}

DE_by_scran <- function(ref_data, data, groups) {
  oracle_calls <- call_DA_scran(ref_data, groups)
  calls <- call_DA_scran(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

calc_DE_discrepancy <- function(ref_data, data, groups, method = "NB") {
  if(method == "scran") {
    DE_calls <- DE_by_scran(ref_data, data, groups)
    oracle_calls <- DE_calls$oracle_calls
    calls <- DE_calls$calls
  } else {
    DE_calls <- DE_by_NB(ref_data, data, groups)
    oracle_calls <- DE_calls$oracle_calls
    calls <- DE_calls$calls
  }

  de <- calls < 0.05
  sim_de <- oracle_calls < 0.05
  
  TP <- sum(de & sim_de)
  FP <- sum(de & !sim_de)
  
  TN <- sum(!de & !sim_de)
  FN <- sum(!de & sim_de)

  return(list(tpr = TP/(TP+FN), fpr = FP/(FP+TN)))  
}

# ------------------------------------------------------------------------------------------------------------
#   Simulation parameters
# ------------------------------------------------------------------------------------------------------------

n <- 5 # replicate number
# Note: 10 seems to be about a minimum within-condition cell/sample number for
# scran marker gene identification
p <- 10000
palette <- generate_highcontrast_palette(p)

ref_data <- "simulated"
# ref_data <- "Morton"
# ref_data <- "Barlow"
# ref_data <- "Athanasiadou_ciona"
# ref_data <- "Athanasiadou_yeast"

k <- 1
asymmetry <- 0.5
sequencing_depth <- 1e5
proportion_da <- 2/3
spike_in <- TRUE
possible_fold_changes <- NULL

iterations <- 200

# ------------------------------------------------------------------------------------------------------------
#   Single-simulation visualizations
# ------------------------------------------------------------------------------------------------------------

if(iterations == 1) {
  
  #   Absolute abundances
  # ------------------------------------------------------------------------------------------------------------
  
  data <- pivot_longer(cbind(sample = 1:nrow(sim_data$abundances), as.data.frame(sim_data$abundances[,1:p])),
                       cols = !sample,
                       names_to = "OTU",
                       values_to = "abundance")
  
  # plot_stacked_bars(data, palette)
  plot_stacked_bars(data, palette, save_name = paste0("01_absolute_",ref_data,".png"))
  
  #   Observed abundances
  # ------------------------------------------------------------------------------------------------------------
  
  # Scenario 1: No correlation between true abundances and observed abundances
  data <- pivot_longer(cbind(sample = 1:nrow(sim_data$observed_counts1), as.data.frame(sim_data$observed_counts1[,1:p])),
                       cols = !sample,
                       names_to = "OTU",
                       values_to = "abundance")
  
  # plot_stacked_bars(data, palette)
  plot_stacked_bars(data, palette, save_name = paste0("02_observed_",ref_data,".png"))
  
  # Scenario 2: Partial correlation between true abundances and observed abundances
  data <- pivot_longer(cbind(sample = 1:nrow(sim_data$observed_counts2), as.data.frame(sim_data$observed_counts2[,1:p])),
                       cols = !sample,
                       names_to = "OTU",
                       values_to = "abundance")
  
  # plot_stacked_bars(data, palette)
  plot_stacked_bars(data, palette, save_name = paste0("02_observed_",ref_data,"_partialcorr.png"))
  
  # Scenario 3: Spike-in normalized
  if(spike_in) {
    spike_in_normalized <- sim_data$observed_counts1[,1:p] / sim_data$observed_counts1[,(p+1)] # row-wise normalize
    data <- pivot_longer(cbind(sample = 1:nrow(sim_data$observed_counts1), as.data.frame(spike_in_normalized)),
                         cols = !sample,
                         names_to = "OTU",
                         values_to = "abundance")
    
    # plot_stacked_bars(data, palette)
    plot_stacked_bars(data, palette, save_name = paste0("02_observed_",ref_data,"_spikedin.png"))
  }
  
  #   Relative abundances
  # ------------------------------------------------------------------------------------------------------------
  
  props <- t(apply(sim_data$abundances, 1, function(x) x/sum(x)))
  data <- pivot_longer(cbind(sample = 1:nrow(props), as.data.frame(props)),
                       cols = !sample,
                       names_to = "OTU",
                       values_to = "abundance") # really, relative abundance
  
  # plot_stacked_bars(data, palette)
  plot_stacked_bars(data, palette, save_name = paste0("03_relative_",ref_data,".png"))

}

# ------------------------------------------------------------------------------------------------------------
#   Visualize fold change distribution across lots of simulations
# ------------------------------------------------------------------------------------------------------------

if(FALSE) {
  fcs <- numeric(iterations)
  for(i in 1:iterations) {
    if(i %% 100 == 0) {
      cat("Iteration",i,"\n")
    }
    sim_data <- simulate_sequence_counts(n = n, p = p, k = k, ref_data = ref_data, asymmetry = asymmetry,
                                         sequencing_depth = sequencing_depth, proportion_da = proportion_da,
                                         spike_in = spike_in, possible_fold_changes = possible_fold_changes)
    totals <- rowSums(sim_data$abundances)
    t1 <- mean(totals[1:n])
    t2 <- mean(totals[(n+1):(n*2)])
    fcs[i] <- t2 / t1
  }
  
  data <- data.frame(fc = fcs)
  ggplot(data, aes(x = fc)) +
    geom_histogram(color = "white") +
    xlim(0, 5)
  ggsave(file.path("images", paste0("00_hist_",ref_data,"_p",p,".png")),
         units = "in",
         dpi = 100,
         height = 5,
         width = 6)
}

# ------------------------------------------------------------------------------------------------------------
#   Run lots of simulations, calculate FPR
# ------------------------------------------------------------------------------------------------------------

plot_data <- data.frame(delta_mean = c(),
                        rate = c(),
                        rate_type = c(),
                        corr = c())

for(i in 1:iterations) {
  if(ref_data == "simulated") {
    # Simulate fresh
    build_simulated_reference(p = p, log_var = 2, log_noise_var = 2)
  }
  sim_data <- simulate_sequence_counts(n = n, p = p, k = k, ref_data = ref_data, asymmetry = asymmetry,
                                       sequencing_depth = sequencing_depth, proportion_da = proportion_da,
                                       spike_in = spike_in, possible_fold_changes = possible_fold_changes)

  # Call differential expression/abundance
  
  # if(i %% 10 == 0) {
  #   cat("Iteration",i,"\n")
  # }
  
  m1 <- mean(rowSums(sim_data$abundances[1:n,]))
  m2 <- mean(rowSums(sim_data$abundances[(n+1):(n*2),]))
  delta_mean_v1 <- max(c(m1, m2)) - min(c(m1, m2))
  delta_mean_v2 <- max(c(m1, m2)) / min(c(m1, m2))
  
  # Discrepancy: true vs. observed
  rates1 <- calc_DE_discrepancy(sim_data$abundances[,1:p], sim_data$observed_counts1[,1:p], sim_data$groups)
  # Discrepancy: true vs. scran-normalized observed
  # rates1a <- calc_DE_discrepancy(sim_data$abundances[,1:p], sim_data$observed_counts1[,1:p], sim_data$groups, method = "scran")
  # Discrepancy: true vs. observed w/ partially informative library sizes
  rates2 <- calc_DE_discrepancy(sim_data$abundances[,1:p], sim_data$observed_counts2[,1:p], sim_data$groups)
  
  cat(paste0("Iteration ",i," FPR: ",round(rates1$fpr, 2),"\n"))
  
  plot_data <- rbind(plot_data,
                     data.frame(delta_mean_v1 = rep(delta_mean_v1, 4), # absolute difference
                                delta_mean_v2 = rep(delta_mean_v2, 4), # fold change
                                rate = c(rates1$fpr, rates1$tpr, rates2$fpr, rates2$tpr),
                                rate_type = rep(c("fpr", "tpr"), 2),
                                corr = c(rep("zero", 2), rep("partial", 2))))
  
}

plot_data$rate_type <- as.factor(plot_data$rate_type)
plot_data$corr <- as.factor(plot_data$corr)

saveRDS(plot_data, file = file.path("output", paste0("simresults_p",p,"_",ref_data,".rds")))

ggplot(plot_data[plot_data$rate_type == "fpr",], aes(x = abs(delta_mean_v2), y = rate, color = corr)) +
  geom_point(size = 2) +
  # geom_smooth(method = "loess") +
  xlim(1, 10) +
  xlab("difference in means")
ggsave(file.path("output", "images", paste0("simresults_p",p,"_",ref_data,".png")),
       units = "in",
       dpi = 100,
       height = 6,
       width = 8)
