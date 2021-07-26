source("path_fix.R")

library(codaDE)
library(uuid)
library(optparse)
library(readr)

# e.g. input = input_100.txt
#      output = res_100 (we'll append row numbers to this)

option_list = list(
  make_option(c("--input"), type = "character", default = NULL,
              help = "input filename", metavar = "character"),
  make_option(c("--output"), type = "character", default = NULL,
              help = "output filename", metavar = "character"),
  make_option(c("--start"), type = "numeric", default = 1,
              help = "start row index", metavar = "numeric"),
  make_option(c("--end"), type = "numeric", default = Inf,
              help = "end row index", metavar = "numeric")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# ------------------------------------------------------------------------------
#   Validate input
# ------------------------------------------------------------------------------

input <- opt$input
output <- opt$output
start <- opt$start
end <- opt$end

if(!file.exists(input)) {
  stop("Input file not found!")
}
if(start < 1 | start > end) {
  stop("Invalid start row!")
}
if(end < 1) {
  stop("Invalid end row!")
}

n <- 10

# These are parameters associated with simulated correlation between features
rhos <- c(0.15, 0.3, 0.45, 0.6)
concentrations <- c(50, 40, 30, 20)
output_fn <- file.path("temp", paste0(output, "_", start, "-", end, ".txt"))

if(!file.exists("temp")) {
  suppressWarnings(dir.create("temp"))
}

# Initialize output
output_file <- file(output_fn)
writeLines(paste0(c("UUID",
                    "P",
                    "CORRP",
                    "LOG_MEAN",
                    "PERTURBATION",
                    "REP_NOISE",
                    "FC_ABSOLUTE",
                    "FC_RELATIVE",
                    "FC_PARTIAL",
                    "BASELINE_CALLS",
                    "MED_ABS_TOTAL",
                    "MED_REL_TOTAL",
                    "PERCENT_DIFF_SIM",
                    "PERCENT_DIFF_REALIZ"), collapse = "\t"),
           output_file)
close(output_file)

# ------------------------------------------------------------------------------
#   Functions
# ------------------------------------------------------------------------------

calc_fc <- function(M) {
  counts_A <- M[1:(nrow(M)/2),]
  counts_B <- M[(nrow(M)/2+1):nrow(M),]
  m1 <- mean(rowSums(counts_A))
  m2 <- mean(rowSums(counts_B))
  m2 / m1
  # max(c(m1, m2)) / min(c(m1, m2))
}

# ------------------------------------------------------------------------------
#   Generate data sets
# ------------------------------------------------------------------------------

# Read chunk to make wishlist
wishlist <- read.table(input)
start <- min(start, nrow(wishlist))
end <- min(end, nrow(wishlist))
wishlist <- wishlist[start:end,]

for(i in 1:nrow(wishlist)) {
  cat(paste0("Processing job ", i, " / ", nrow(wishlist), "\n"))
  
  # Parse data set
  job <- wishlist[i,]
  
  half_p <- round(job$P/2) # need this repeatedly for later calculations
  
  base_correlation <- matrix(0, job$P, job$P)
  if(job$CORRP == 0) {
    # Independent features
    base_correlation <- diag(job$P)
    concentration <- 1e6
  } else {
    # Partial correlation
    base_correlation[1:half_p,1:half_p] <- rhos[job$CORRP]
    diag(base_correlation) <- 1
    concentration <- job$P + concentrations[job$CORRP]
  }
  
  # Create reference distributions
  data_obj <- build_simulated_reference(p = job$P,
                                        log_mean = job$LOG_MEAN,
                                        perturb_size = job$PERTURBATION,
                                        base_correlation = base_correlation,
                                        concentration = concentration)
  
  # Sample from these to create the data set
  sim_data <- simulate_sequence_counts(n = n,
                                       p = job$P,
                                       data_obj = data_obj,
                                       replicate_noise = job$REP_NOISE,
                                       proportion_da = job$PERCENT_DIFF)
  
  # Visualize via
  # plot_stacked_bars(sim_data$abundances)

  # Generate an ID for this data set
  uuid <- UUIDgenerate()
  
  fc_abs <- calc_fc(sim_data$abundances)
  fc_rel <- calc_fc(sim_data$observed_counts1)
  fc_par <- calc_fc(sim_data$observed_counts2)
  
  calls <- call_DA_NB(sim_data$abundances, sim_data$groups)$pval
  
  med_abs <- mean(rowSums(sim_data$abundances))
  med_rel <- mean(rowSums(sim_data$observed_counts1))
  
  percent_diff_realiz <- sum(p.adjust(calls, method = "BH") < 0.05) / length(calls)
  
  # baseline_calls <- ifelse(calls < 0.05, 0, 1)
  # baseline_calls_mtc <- ifelse(calls < 0.05/p, 0, 1)
  
  # Save it to "datasets" folder
  output_dir <- file.path("output", "datasets")
  dir.create(output_dir, showWarnings = FALSE)
  saveRDS(list(reference = data_obj,
               simulation = sim_data),
          file = file.path(output_dir, paste0(uuid, ".rds")))
  
  # Add to output file
  write_delim(data.frame(UUID = uuid,
                         P = job$P,
                         CORRP = job$CORRP,
                         LOG_MEAN = job$LOG_MEAN,
                         PERTURBATION = job$PERTURBATION,
                         REP_NOISE = job$REP_NOISE,
                         FC_ABSOLUTE = fc_abs,
                         FC_RELATIVE = fc_rel,
                         FC_PARTIAL = fc_par,
                         BASELINE_CALLS = paste0(round(calls, 10), collapse = ";"),
                         MED_ABS_TOTAL = med_abs,
                         MED_REL_TOTAL = med_rel,
                         PERCENT_DIFF_SIM = job$PERCENT_DIFF,
                         PERCENT_DIFF_REALIZ = percent_diff_realiz),
              output_fn,
              delim = "\t",
              append = TRUE)
}

