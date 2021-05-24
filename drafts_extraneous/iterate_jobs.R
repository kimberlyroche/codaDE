source("path_fix.R")

library(tidyverse)
library(gridExtra)
library(codaDE)
library(uuid)
library(optparse)

option_list = list(
  make_option(c("--p"), type = "numeric", default = 100,
              help = "number of genes", metavar = "numeric"),
  make_option(c("--corrp"), type = "numeric", default = 0,
              help = "proportion of features to simulate as moderately positively correlated", metavar = "numeric"),
  make_option(c("--i"), type = "numeric", default = 20,
              help = "number of iterations", metavar = "numeric")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# ------------------------------------------------------------------------------------------------------------
#   Simulation parameters: CHANGE THESE
# ------------------------------------------------------------------------------------------------------------

# This stuff changes!
p <- opt$p
correlation_prop <- opt$corrp
# This is the proportion of features correlated: 0, 0.5, 1
iterations <- opt$i

# ------------------------------------------------------------------------------------------------------------
#   Functions
# ------------------------------------------------------------------------------------------------------------

DE_by_NB <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_NB(ref_data, groups)
  }
  calls <- call_DA_NB(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

DE_by_DESeq2 <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_DESeq2(ref_data, groups)
  }
  calls <- call_DA_DESeq2(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

DE_by_MAST <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_MAST(ref_data, groups)
  }
  calls <- call_DA_MAST(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

DE_by_ALDEx2 <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_ALDEx2(ref_data, groups)
  }
  calls <- call_DA_ALDEx2(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

DE_by_scran <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_scran(ref_data, groups)
  }
  calls <- call_DA_scran(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

calc_DE_discrepancy <- function(ref_data, data, groups, method = "NB",
                                oracle_calls = NULL) {
  if(method == "NB") {
    DE_calls <- DE_by_NB(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
  if(method == "DESeq2") {
    DE_calls <- DE_by_DESeq2(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
  if(method == "MAST") {
    DE_calls <- DE_by_MAST(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
  if(method == "ALDEx2") {
    DE_calls <- DE_by_ALDEx2(ref_data, data, groups,
                             oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
  if(method == "scran") {
    DE_calls <- DE_by_scran(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
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
#   Simulation parameters: CONSTANT
# ------------------------------------------------------------------------------------------------------------

output_dir <- paste0("p", p, "_corrp", correlation_prop)
dir.create(file.path("output", output_dir), showWarnings = FALSE)

if(correlation_prop == 0) {
  base_correlation <- diag(p)
  concentration <- 1e6
} else if(correlation_prop == 1) {
  base_correlation <- matrix(0.5, p, p)
  diag(base_correlation) <- 1
  concentration <- p + 100
} else {
  # Assume 50%
  half_p <- round(p/2)
  base_correlation <- matrix(0, p, p)
  base_correlation[1:half_p,1:half_p] <- 0.5
  diag(base_correlation) <- 1
  concentration <- p + 100
}

n <- 10 # replicate number
asymmetry <- 1
proportion_da <- 0.75
spike_in <- FALSE

# ------------------------------------------------------------------------------------------------------------
#   Run lots of simulations, calculate FPR
# ------------------------------------------------------------------------------------------------------------

plot_data <- data.frame(uuid = c(),
                        delta_mean_v1 = c(),
                        delta_mean_v2 = c(),
                        # cor_totals = c(),
                        gap_totals = c(),
                        rate = c(),
                        rate_type = c(),
                        method = c(),
                        partial_info = c())

for(i in 1:iterations) {

  cat(paste0("------------ STARTING ITERATION #", i, "\n"))

  uuid <- UUIDgenerate()
  data_obj <- build_simulated_reference(p = p,
                                        log_mean = 1,
                                        log_var = 2,
                                        log_noise_var = 2,
                                        base_correlation = base_correlation,
                                        concentration = concentration)
  sim_data <- simulate_sequence_counts(n = n,
                                       p = p,
                                       data_obj = data_obj,
                                       asymmetry = asymmetry,
                                       proportion_da = proportion_da,
                                       spike_in = spike_in)
  saveRDS(sim_data, file = file.path("output", output_dir, paste0(uuid, ".rds")))

  # Characterize the data
  
  r1 <- rowSums(sim_data$abundances[1:n,])
  r2 <- rowSums(sim_data$abundances[(n+1):(n*2),])
  m1 <- mean(r1)
  m2 <- mean(r2)
  absolute_diff <- max(c(m1, m2)) - min(c(m1, m2)) # absolute difference
  fold_diff <- max(c(m1, m2)) / min(c(m1, m2)) # absolute fold change
  
  # Get NB calls first, as a baseline ("oracle")
  cat(paste0("\tEvaluating method: NB\n"))
  rates_baseline <- calc_DE_discrepancy(sim_data$abundances[,1:p],
                                        sim_data$observed_counts1[,1:p],
                                        sim_data$groups)
  oracle_calls <- rates_baseline$oracle_calls
  plot_data <- rbind(plot_data,
                     data.frame(uuid = uuid,
                                absolute_diff = rep(absolute_diff, 2),
                                fold_diff = rep(fold_diff, 2),
                                rate = c(rates_baseline$fpr, rates_baseline$tpr),
                                rate_type = c("fpr", "tpr"),
                                method = rep("baseline", 2),
                                partial_info = FALSE))
  
  # Get NB calls on partially informative total abundances
  rates_partial <- calc_DE_discrepancy(sim_data$abundances[,1:p],
                                       sim_data$observed_counts2[,1:p],
                                       sim_data$groups)
  plot_data <- rbind(plot_data,
                     data.frame(uuid = uuid,
                                absolute_diff = rep(absolute_diff, 2),
                                fold_diff = rep(fold_diff, 2),
                                rate = c(rates_partial$fpr, rates_partial$tpr),
                                rate_type = c("fpr", "tpr"),
                                method = rep("baseline", 2),
                                partial_info = TRUE))

  methods <- c("DESeq2", "scran", "MAST", "ALDEx2")
  other_methods_for_partials <- c("DESeq2")
  for(method in methods) {
    cat(paste0("\tEvaluating method: ", method, "\n"))
    rates <- calc_DE_discrepancy(sim_data$abundances[,1:p],
                                 sim_data$observed_counts1[,1:p],
                                 sim_data$groups,
                                 method = method,
                                 oracle_calls = oracle_calls)
    
    plot_data <- rbind(plot_data,
                       data.frame(uuid = uuid,
                                  absolute_diff = rep(absolute_diff, 2),
                                  fold_diff = rep(fold_diff, 2),
                                  rate = c(rates$fpr, rates$tpr),
                                  rate_type = c("fpr", "tpr"),
                                  method = c(rep(method, 2)),
                                  partial_info = FALSE))
    
    if(method %in% other_methods_for_partials) {
      rates_partial <- calc_DE_discrepancy(sim_data$abundances[,1:p],
                                   sim_data$observed_counts2[,1:p],
                                   sim_data$groups,
                                   method = method,
                                   oracle_calls = oracle_calls)
      
      plot_data <- rbind(plot_data,
                         data.frame(uuid = uuid,
                                    absolute_diff = rep(absolute_diff, 2),
                                    fold_diff = rep(fold_diff, 2),
                                    rate = c(rates_partial$fpr, rates_partial$tpr),
                                    rate_type = c("fpr", "tpr"),
                                    method = c(rep(method, 2)),
                                    partial_info = TRUE))
    }
    
  }
}

plot_data$rate_type <- as.factor(plot_data$rate_type)
plot_data$method <- as.factor(plot_data$method)

saveRDS(plot_data, file = file.path("output", output_dir, paste0("simresults_p",p,"_simulated_",UUIDgenerate(),".rds")))
