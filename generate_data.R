source("path_fix.R")

library(RSQLite)
library(tidyverse)
library(codaDE)
library(uuid)
library(optparse)

option_list = list(
  make_option(c("--p"), type = "numeric", default = 100,
              help = "number of genes", metavar = "numeric"),
  make_option(c("--corrp"), type = "numeric", default = 0,
              help = "simulate net positively correlated features", metavar = "numeric"),
  make_option(c("--iter_start"), type = "numeric", default = 1,
              help = "index of first iteration", metavar = "numeric"),
  make_option(c("--iter_end"), type = "numeric", default = 100,
              help = "index of last iteration", metavar = "numeric")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

p <- opt$p
corrp <- opt$corrp
iter_start <- opt$iter_start
iter_end <- opt$iter_end

if(iter_start < 1 | iter_start > iter_end | iter_end > 100) {
  stop("Iteration start and end must be in range 1..100!\n")
}

calc_fc <- function(M) {
  counts_A <- M[1:(nrow(M)/2),]
  counts_B <- M[(nrow(M)/2+1):nrow(M),]
  m1 <- mean(rowSums(counts_A))
  m2 <- mean(rowSums(counts_B))
  max(c(m1, m2)) / min(c(m1, m2))
}

calc_fc <- function(M) {
  counts_A <- M[1:(nrow(M)/2),]
  counts_B <- M[(nrow(M)/2+1):nrow(M),]
  m1 <- mean(rowSums(counts_A))
  m2 <- mean(rowSums(counts_B))
  max(c(m1, m2)) / min(c(m1, m2))
}

# ------------------------------------------------------------------------------
#   Simulate (and save) data set
# ------------------------------------------------------------------------------

# Set up simulation parameters
if(corrp == 1) {
  # Generate roughly 50% positively correlated features. Remaining features will
  # be random (and roughly symmetrical) in their +/- correlation.
  half_p <- round(p/2)
  base_correlation <- matrix(0, p, p)
  base_correlation[1:half_p,1:half_p] <- 0.8
  diag(base_correlation) <- 1
  concentration <- p + 10
} else {
  base_correlation <- diag(p)
  concentration <- 1e6
}

n <- 10 # replicate number
asymmetry <- 1
proportion_da <- 0.75
spike_in <- FALSE

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
# Increase the "busy" timeout; default is too short
discard <- dbExecute(conn, "PRAGMA busy_timeout = 60000;")

perturbations <- seq(from = 0.1, to = 4, length.out = 100)
for(i in iter_start:iter_end) {
  cat(paste0("Iteration ", i, " / ", iter_end, "\n"))
  uuid <- UUIDgenerate()
  
  # Create data set
  data_obj <- build_simulated_reference(p = p,
                                        log_mean = 1,
                                        log_noise_var = perturbations[i],
                                        base_correlation = base_correlation,
                                        concentration = concentration)
  sim_data <- simulate_sequence_counts(n = n,
                                       p = p,
                                       data_obj = data_obj,
                                       asymmetry = asymmetry,
                                       proportion_da = proportion_da,
                                       spike_in = spike_in)

  prop_de <- sum(calc_threshold_DA(sim_data$abundances) == 0) / ncol(sim_data$abundances)
  fc_abs <- calc_fc(sim_data$abundances)
  fc_rel <- calc_fc(sim_data$observed_counts1)
  fc_par <- calc_fc(sim_data$observed_counts2)
  
  # Save it to "datasets" folder
  output_dir <- file.path("output", "datasets")
  dir.create(output_dir, showWarnings = FALSE)
  saveRDS(list(reference = data_obj,
               simulation = sim_data),
          file = file.path(output_dir, paste0(uuid, ".rds")))
  
  # ------------------------------------------------------------------------------
  #   Add data set to DB
  # ------------------------------------------------------------------------------
  res <- dbExecute(conn, paste0("INSERT into datasets(UUID,P,CORRP,FC_ABSOLUTE,FC_RELATIVE,FC_PARTIAL,PERCENT_DIFF,TIMESTAMP) ",
                                "VALUES(",
                                "'",uuid,"',",
                                p,",",
                                as.integer(corrp),",",
                                fc_abs,",",
                                fc_rel,",",
                                fc_par,",",
                                prop_de,",",
                                "'",get_timestamp(),"'",
                                ")"))
}
# dbGetQuery(conn, "SELECT * FROM datasets")
dbDisconnect(conn)
