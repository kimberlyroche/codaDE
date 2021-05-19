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
              help = "proportion of features to simulate as moderately positively correlated", metavar = "numeric"),
  make_option(c("--iter"), type = "numeric", default = 1,
              help = "number of data sets to generate", metavar = "numeric")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

p <- opt$p
corrp <- opt$corrp
iter <- opt$iter

# ------------------------------------------------------------------------------
#   Simulate (and save) data set
# ------------------------------------------------------------------------------

# Set up simulation parameters
if(corrp == 0) {
  base_correlation <- diag(p)
  concentration <- 1e6
} else if(corrp == 1) {
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

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
for(i in 1:iter) {
  uuid <- UUIDgenerate()
  
  # Create data set
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
  
  # Save it to "datasets" folder
  output_dir <- file.path("output", "datasets")
  dir.create(output_dir, showWarnings = FALSE)
  saveRDS(list(reference = data_obj,
               simulation = sim_data),
          file = file.path(output_dir, paste0(uuid, ".rds")))
  
  # ------------------------------------------------------------------------------
  #   Add data set to DB
  # ------------------------------------------------------------------------------
  res <- dbExecute(conn, paste0("INSERT into datasets(UUID,P,CORRP,TIMESTAMP) ",
                                "VALUES(",
                                "'",uuid,"',",
                                p,",",
                                corrp,",",
                                "'",get_timestamp(),"'",
                                ")"))
}
# dbGetQuery(conn, "SELECT * FROM datasets")
dbDisconnect(conn)


