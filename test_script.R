# 2021-02-27
# This is a constantly-edited test script for generating simulated data sets and evaluating their properties.

if(R.version$major == 4 & Sys.info()[[1]] != "Windows") {
  cat("Updating lib paths...\n")
  .libPaths(c("/gpfs/fs1/data/mukherjeelab/roche/Rlibs", .libPaths()[2]))
}

library(tidyverse)
library(gridExtra)
library(codaDE)
library(uuid)

# ------------------------------------------------------------------------------------------------------------
#   Functions
# ------------------------------------------------------------------------------------------------------------

DE_by_NB <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- sapply(1:p, function(idx) { call_DA_NB(ref_data, groups, idx) } )
  }
  calls <- sapply(1:p, function(idx) { call_DA_NB(data, groups, idx) } )
  return(list(oracle_calls = oracle_calls, calls = calls))
}

DE_by_edgeR <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_edgeR(ref_data, groups)
  }
  calls <- call_DA_edgeR(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

DE_by_Seurat <- function(ref_data, data, groups, method, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_Seurat(ref_data, groups, method = method)
  }
  calls <- call_DA_Seurat(data, groups, method = method)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

DE_by_MAST <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_MAST(ref_data, groups)
  }
  calls <- call_DA_MAST(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

DE_by_scran <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_scran(ref_data, groups)
  }
  calls <- call_DA_scran(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

DE_by_ALDEx2 <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_ALDEx2(ref_data, groups)
  }
  calls <- call_DA_ALDEx2(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

calc_DE_discrepancy <- function(ref_data, data, groups, method = "NB",
                                oracle_calls = NULL) {
  if(method == "scran") {
    DE_calls <- DE_by_scran(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls
    }
    calls <- DE_calls$calls
  } else if(method == "edgeR") {
    DE_calls <- DE_by_edgeR(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls
    }
    calls <- DE_calls$calls
  } else if(method == "wilcox" | method == "DESeq2") {
    DE_calls <- DE_by_Seurat(ref_data, data, groups, method = method, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls
    }
    calls <- DE_calls$calls
  } else if(method == "MAST") {
    DE_calls <- DE_by_MAST(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls
    }
    calls <- DE_calls$calls
  } else if(method == "ALDEx2") {
    DE_calls <- DE_by_ALDEx2(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls
    }
    calls <- DE_calls$calls
  } else {
    DE_calls <- DE_by_NB(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls
    }
    calls <- DE_calls$calls
  }

  de <- calls < 0.05
  sim_de <- oracle_calls < 0.05
  
  TP <- sum(de & sim_de)
  FP <- sum(de & !sim_de)
  
  TN <- sum(!de & !sim_de)
  FN <- sum(!de & sim_de)

  return(list(tpr = TP/(TP+FN), fpr = FP/(FP+TN), oracle_calls = oracle_calls))  
}

# ------------------------------------------------------------------------------------------------------------
#   Simulation parameters
# ------------------------------------------------------------------------------------------------------------

n <- 10 # replicate number
# Note: 10 seems to be about a minimum within-condition cell/sample number for
# scran marker gene identification

# ref_data <- "simulated_bulk"
# ref_data <- "simulated_16S"
# ref_data <- "simulated_sc"
# ref_data <- "Morton"
# ref_data <- "Barlow"
# ref_data <- "Athanasiadou_ciona"
# ref_data <- "Athanasiadou_yeast"

p <- 1000

if(!exists("p")) {
  if(ref_data %in% c("simulated_bulk", "Athanasiadou_ciona", "Athanasiadou_yeast")) {
    p <- 15000
  } else if(ref_data %in% c("simulated_16S", "Morton", "Barlow")) {
    p <- 1000
  } else if(ref_data %in% c("simulated_sc")) {
    p <- 5000
  }
}

palette <- generate_highcontrast_palette(p)

asymmetry <- 1
proportion_da <- 0.75
spike_in <- TRUE

iterations <- 20

# ------------------------------------------------------------------------------------------------------------
#   Single-simulation visualizations
# ------------------------------------------------------------------------------------------------------------

if(iterations == 1) {
  
  # Simulate fresh
  data_obj <- build_simulated_reference(p = p, log_mean = 1, log_var = 2, log_noise_var = 2)
  sim_data <- simulate_sequence_counts(n = n,
                                       p = p,
                                       data_obj = data_obj,
                                       asymmetry = asymmetry,
                                       proportion_da = proportion_da,
                                       spike_in = spike_in)
  
  m1 <- mean(rowSums(sim_data$abundances[1:n,]))
  m2 <- mean(rowSums(sim_data$abundances[(n+1):(n*2),]))
  fc <- max(c(m1, m2)) / min(c(m1, m2))
  
  cat("Fold change:",round(fc, 2),"\n")
  
  #   Absolute abundances
  # ------------------------------------------------------------------------------------------------------------
  
  data <- pivot_longer(cbind(sample = 1:nrow(sim_data$abundances), as.data.frame(sim_data$abundances[,1:p])),
                       cols = !sample,
                       names_to = "feature",
                       values_to = "abundance")

  # plot_stacked_bars(data, palette, save_name = "inset_3.png")
    
  plot_stacked_bars(data, palette)
  plot_stacked_bars(data, palette, save_name = paste0("01_absolute_",ref_data,".png"))
  
  #   Observed abundances
  # ------------------------------------------------------------------------------------------------------------
  
  # Scenario 1: No correlation between true abundances and observed abundances
  data <- pivot_longer(cbind(sample = 1:nrow(sim_data$observed_counts1), as.data.frame(sim_data$observed_counts1[,1:p])),
                       cols = !sample,
                       names_to = "feature",
                       values_to = "abundance")
  
  plot_stacked_bars(data, palette)
  plot_stacked_bars(data, palette, save_name = paste0("02_observed_",ref_data,".png"))
  
  # Scenario 2: Partial correlation between true abundances and observed abundances
  data <- pivot_longer(cbind(sample = 1:nrow(sim_data$observed_counts2), as.data.frame(sim_data$observed_counts2[,1:p])),
                       cols = !sample,
                       names_to = "feature",
                       values_to = "abundance")
  
  plot_stacked_bars(data, palette)
  plot_stacked_bars(data, palette, save_name = paste0("02_observed_",ref_data,"_partialcorr.png"))
  
  # Scenario 3: Spike-in normalized
  if(spike_in) {
    spike_in_normalized <- sim_data$observed_counts1[,1:p] / sim_data$observed_counts1[,(p+1)] # row-wise normalize
    data <- pivot_longer(cbind(sample = 1:nrow(sim_data$observed_counts1), as.data.frame(spike_in_normalized)),
                         cols = !sample,
                         names_to = "feature",
                         values_to = "abundance")
    
    plot_stacked_bars(data, palette)
    plot_stacked_bars(data, palette, save_name = paste0("02_observed_",ref_data,"_spikedin.png"))
  }
  
  #   Relative abundances
  # ------------------------------------------------------------------------------------------------------------
  
  props <- t(apply(sim_data$abundances, 1, function(x) x/sum(x)))
  data <- pivot_longer(cbind(sample = 1:nrow(props), as.data.frame(props)),
                       cols = !sample,
                       names_to = "feature",
                       values_to = "abundance") # really, relative abundance
  
  plot_stacked_bars(data, palette)
  plot_stacked_bars(data, palette, save_name = paste0("03_relative_",ref_data,".png"))

}

# ------------------------------------------------------------------------------------------------------------
#   Visualize fold change distribution across lots of simulations
# ------------------------------------------------------------------------------------------------------------

if(FALSE) {
  fcs <- numeric(iterations)
  data_obj <- build_simulated_reference(p = p, log_mean = 1, log_var = 2, log_noise_var = 2)
  for(i in 1:iterations) {
    if(i %% 100 == 0) {
      cat("Iteration",i,"\n")
    }
    # Resample the same data
    sim_data <- simulate_sequence_counts(n = n,
                                         p = p,
                                         data_obj = data_obj,
                                         asymmetry = asymmetry,
                                         proportion_da = proportion_da,
                                         spike_in = spike_in)
    totals <- rowSums(sim_data$abundances)
    t1 <- mean(totals[1:n])
    t2 <- mean(totals[(n+1):(n*2)])
    fcs[i] <- t2 / t1
  }
  
  data <- data.frame(fc = fcs)
  ggplot(data, aes(x = fc)) +
    geom_histogram(color = "white") +
    xlim(0, 5)
  ggsave(file.path("output", "images", paste0("00_hist_",ref_data,"_p",p,".png")),
         units = "in",
         dpi = 100,
         height = 5,
         width = 6)
}

# ------------------------------------------------------------------------------------------------------------
#   Run lots of simulations, calculate FPR
# ------------------------------------------------------------------------------------------------------------

plot_data <- data.frame(delta_mean_v1 = c(),
                        delta_mean_v2 = c(),
                        cor_totals = c(),
                        rate = c(),
                        rate_type = c(),
                        method = c())

for(i in 1:iterations) {
  cat("------------ STARTING ITERATION",i,"\n")
  # if(ref_data == "simulated_bulk") {
  #   build_simulated_reference(p = p, log_mean = 0, log_var = 2, log_noise_var = 2, save_name = ref_data)
  # } else if(ref_data == "simulated_16S") {
  #   build_simulated_reference(p = p, log_mean = -1, log_var = 3, log_noise_var = 2, save_name = ref_data)
  # } else if(ref_data == "simulated_sc") {
  #   build_simulated_reference(p = p, log_mean = -1, log_var = 2, log_noise_var = 1, save_name = ref_data)
  # }
  data_obj <- build_simulated_reference(p = p, log_mean = 1, log_var = 2, log_noise_var = 2)
  sim_data <- simulate_sequence_counts(n = n,
                                       p = p,
                                       data_obj = data_obj,
                                       asymmetry = asymmetry,
                                       proportion_da = proportion_da,
                                       spike_in = spike_in)

  # Call differential expression/abundance
  
  # if(i %% 10 == 0) {
  #   cat("Iteration",i,"\n")
  # }

  r1 <- rowSums(sim_data$abundances[1:n,])
  r2 <- rowSums(sim_data$abundances[(n+1):(n*2),])
  m1 <- mean(r1)
  m2 <- mean(r2)
  delta_mean_v1 <- max(c(m1, m2)) - min(c(m1, m2))
  delta_mean_v2 <- max(c(m1, m2)) / min(c(m1, m2))
  
  totals <- rowSums(sim_data$abundances)
  totals1 <- rowSums(sim_data$observed_counts1)
  totals2 <- rowSums(sim_data$observed_counts2)
  
  # Discrepancy: true vs. observed
  # NB model first
  rates_baseline <- calc_DE_discrepancy(sim_data$abundances[,1:p], sim_data$observed_counts1[,1:p], sim_data$groups)
  # Save the oracle calls as our gold standard for this simulated data set
  oracle_calls <- rates_baseline$oracle_calls
  plot_data <- rbind(plot_data,
                     data.frame(delta_mean_v1 = rep(delta_mean_v1, 2),
                                delta_mean_v2 = rep(delta_mean_v2, 2),
                                cor_totals = cor(totals, totals1),
                                rate = c(rates_baseline$fpr, rates_baseline$tpr),
                                rate_type = c("fpr", "tpr"),
                                method = rep("baseline", 2)))
  
  rates_partial <- calc_DE_discrepancy(sim_data$abundances[,1:p], sim_data$observed_counts2[,1:p], sim_data$groups)
  plot_data <- rbind(plot_data,
                     data.frame(delta_mean_v1 = rep(delta_mean_v1, 2),
                                delta_mean_v2 = rep(delta_mean_v2, 2),
                                cor_totals = cor(totals, totals2),
                                rate = c(rates_partial$fpr, rates_partial$tpr),
                                rate_type = c("fpr", "tpr"),
                                method = rep("spike_in", 2)))

  methods <- c("DESeq2", "scran", "MAST", "ALDEx2", "edgeR_TMM")
  for(method in methods) {
    rates <- calc_DE_discrepancy(sim_data$abundances[,1:p],
                                 sim_data$observed_counts1[,1:p],
                                 sim_data$groups,
                                 method = method,
                                 oracle_calls = oracle_calls)
    
    plot_data <- rbind(plot_data,
                       data.frame(delta_mean_v1 = rep(delta_mean_v1, 2), # absolute difference
                                  delta_mean_v2 = rep(delta_mean_v2, 2), # fold change
                                  cor_totals = cor(totals, totals1),
                                  rate = c(rates$fpr, rates$tpr),
                                  rate_type = c("fpr", "tpr"),
                                  method = c(rep(method, 2))))
    
    # cat(paste0("Iteration ",i," FPR (",method,"): ",round(rates$fpr, 2),"\n"))
  }
}

plot_data$rate_type <- as.factor(plot_data$rate_type)
plot_data$method <- as.factor(plot_data$method)

saveRDS(plot_data, file = file.path("output", paste0("simresults_p",p,"_",ref_data,"_",UUIDgenerate(),".rds")))

# ggplot(plot_data[plot_data$rate_type == "fpr",], aes(x = abs(delta_mean_v2), y = rate, color = method)) +
#   geom_point(size = 2) +
#   # geom_smooth(method = "loess") +
#   xlim(1, 10) +
#   xlab("difference in means")
# ggsave(file.path("output", "images", paste0("simresults_p",p,"_",ref_data,".png")),
#        units = "in",
#        dpi = 100,
#        height = 6,
#        width = 8)
