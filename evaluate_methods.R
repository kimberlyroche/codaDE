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
results <- NULL

for(i in 1:nrow(wishlist)) {
  cat(paste0("Processing job ", i, " / ", nrow(wishlist), "\n"))

  # Parse data set
  job <- wishlist[i,]
  data <- readRDS(file.path(output_dir, paste0(job$uuid, ".rds")))

  # Get baseline differential abundance calls
  if(job$baseline == "oracle") {
    oracle_calls <- as.numeric(strsplit(baseline_calls %>% filter(UUID == job$uuid) %>% pull(BASELINE_CALLS), ";")[[1]])
  } else if(job$baseline == "self") {
    oracle_calls <- NULL
  }

  if(job$partial_info == 0) {
    rates <- calc_DE_discrepancy(data$simulation$abundances,
                                 data$simulation$observed_counts1,
                                 data$simulation$groups,
                                 method = job$method,
                                 oracle_calls = oracle_calls)
  } else if(job$partial_info == 1) {
    rates <- calc_DE_discrepancy(data$simulation$abundances,
                                 data$simulation$observed_counts2,
                                 data$simulation$groups,
                                 method = job$method,
                                 oracle_calls = oracle_calls)
  } 

  if(is.null(results)) {
    results <- cbind(job, method = method, tpr = rates$tpr, fpr = rates$fpr)
  } else {
    results <- rbind(results,
      cbind(job, method = method, tpr = rates$tpr, fpr = rates$fpr))
  }
}

output_fn <- file.path("temp", paste0(output, "_", start, "-", end, ".txt"))
write.table(results, file = output_fn)
