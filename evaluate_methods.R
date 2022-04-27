source("path_fix.R")

library(tidyverse)
library(codaDE)
library(RSQLite)
library(optparse)

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
#  Validate input
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

output_fn <- file.path("temp", paste0(output, "_", start, "-", end, ".txt"))

# Initialize output
output_file <- file(output_fn)
writeLines(paste0(c("ID",
                    "uuid",
                    "baseline",
                    "partial_info",
                    "method",
                    "observed_type",
                    "baseline_calls",
                    "baseline_betas",
                    "calls",
                    "betas"), collapse = "\t"),
           output_file)
close(output_file)

# ------------------------------------------------------------------------------
#  Pull jobs to run, then run 'em!
# ------------------------------------------------------------------------------

wishlist <- read.table(input)
start <- min(start, nrow(wishlist))
end <- min(end, nrow(wishlist))
wishlist <- wishlist[start:end,]

if(any(!(wishlist$baseline %in% c("self", "oracle")))) {
  stop("Invalid reference for some job(s)!")
}
if(any(!(wishlist$partial_info %in% c(0, 1)))) {
  stop("Invalid partial information flag for some job(s)!")
}

baseline_calls <- NULL
if(any(wishlist$baseline == "oracle")) {
  # Pull baseline calls
  conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
  baseline_calls <- dbGetQuery(conn, paste0("SELECT UUID, BASELINE_CALLS FROM datasets WHERE UUID IN ('",
    paste0(unique(wishlist$uuid), collapse = "', '"),"')"))
  dbDisconnect(conn)
}

output_dir <- file.path("output", "datasets")

for(i in 1:nrow(wishlist)) {
  cat(paste0("Processing job ", i, " / ", nrow(wishlist), "\n"))

  # Parse data set
  job <- wishlist[i,]
  data <- readRDS(file.path(output_dir, paste0(job$uuid, ".rds")))

  if(job$observed_type == "cpm") {
    data$simulation$observed_counts1 <- t(apply(data$simulation$observed_counts1, 1, function(x) x/sum(x)))
    data$simulation$observed_counts1 <- data$simulation$observed_counts1*1e06
    data$simulation$observed_counts1 <- round(data$simulation$observed_counts1)
  }

  # Get baseline differential abundance calls
  if(job$baseline == "oracle") {
    oracle_calls <- as.numeric(strsplit(baseline_calls %>% filter(UUID == job$uuid) %>% pull(BASELINE_CALLS), ";")[[1]])
  } else if(job$baseline == "self") {
    oracle_calls <- NULL
  }

  if(job$method == "DESeq2_control") {
    # Find "housekeeping" gene analogs
    log_ab <- log(t(data$simulation$abundances) + 0.5)
    covar <- apply(log_ab, 1, function(x) sd(x)/mean(x))
    bottom10 <- quantile(abs(covar), probs = c(0.1))
    control_indices <- sample(which(abs(covar) < bottom10), size = min(sum(abs(covar) < bottom10), 10))
    job$method <- "DESeq2"
  } else {
    control_indices <- NULL
  }

  if(job$partial_info == 0) {
    all_calls <- tryCatch({
      DA_wrapper(data$simulation$abundances,
                 data$simulation$observed_counts1,
                 data$simulation$groups,
                 method = job$method,
                 oracle_calls = oracle_calls,
                 control_indices = control_indices)
    },
    error = function(cond) { NULL },
    warning = function(cond) { NULL },
    finally = { })
  } else if(job$partial_info == 1) {
    all_calls <- tryCatch({
      DA_wrapper(data$simulation$abundances,
                 data$simulation$observed_counts2,
                 data$simulation$groups,
                 method = job$method,
                 oracle_calls = oracle_calls,
                 control_indices = control_indices)
    },
    error = function(cond) { NULL },
    warning = function(cond) { NULL },
    finally = { })
  }

  if(!is.null(control_indices)) {
    job$method <- "DESeq2_control"
  }
  
  if(!is.null(all_calls)) {
    if(job$baseline == "oracle") {
      # The baseline calls are already present in the datasets table. In the interest of not
      # keeping multiple copies of these - which could get out of sync - we'll force ourselves
      # to refer back to the datasets table to find the oracle calls.
      results_row <- cbind(job,
                           baseline_calls = NA,
                           baseline_betas = NA,
                           calls = paste0(round(all_calls$calls$pval, 10), collapse = ";"),
                           betas = paste0(round(all_calls$calls$beta, 10), collapse = ";"))
    } else {
      results_row <- cbind(job,
                           baseline_calls = paste0(round(all_calls$oracle_calls$pval, 10), collapse = ";"),
                           baseline_betas = paste0(round(all_calls$oracle_calls$beta, 10), collapse = ";"),
                           calls = paste0(round(all_calls$calls$pval, 10), collapse = ";"),
                           betas = paste0(round(all_calls$calls$beta, 10), collapse = ";"))
    }

    # Add to output file
    write_delim(results_row, output_fn, delim = "\t", append = TRUE)
  }
}
