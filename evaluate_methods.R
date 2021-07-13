source("path_fix.R")

library(tidyverse)
library(codaDE)
library(optparse)

option_list = list(
  make_option(c("--method"), type = "character", default = NULL,
              help = "differential abundance method", metavar = "character"),
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

allowed_methods <- c("ALDEx2", "DESeq2", "MAST", "scran")

# ------------------------------------------------------------------------------
#  Validate input
# ------------------------------------------------------------------------------

method <- opt$method
input <- opt$input
output <- opt$output
start <- opt$start
end <- opt$end

if(!(method %in% allowed_methods)) {
  stop("Unrecognized method!")
}
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

if(any(!(wishlist$baseline %in% c("self", "threshold")))) {
  stop("Invalid reference for some job(s)!")
}
if(any(!(wishlist$partial_info %in% c(0, 1)))) {
  stop("Invalid partial information flag for some job(s)!")
}

output_dir <- file.path("output", "datasets")
results <- NULL

for(i in 1:nrow(wishlist)) {
  cat(paste0("Processing job ", i, " / ", nrow(wishlist), "\n"))

  # Parse data set
  job <- wishlist[i,]
  data <- readRDS(file.path(output_dir, paste0(job$uuid, ".rds")))
  
  # LEFT OFF HERE 7/13; NEED TO PARSE BASELINE CALLS FROM DATASETS DB TABLE
  # FOR USE AS ORACLE!!!

  # Get baseline differential abundance calls
  if(job$baseline == "threshold") {
    # "Differential" features will be those with mean fold change in abundance of 
    # >= 1.5 or <= 0.5
    oracle_calls <- calc_threshold_DA(data$simulation$abundances)
  } else if(job$baseline == "self") {
    oracle_calls <- NULL
  }

  if(job$partial_info == 0) {
    rates <- calc_DE_discrepancy(data$simulation$abundances,
                                 data$simulation$observed_counts1,
                                 data$simulation$groups,
                                 method = method,
                                 oracle_calls = oracle_calls)
  } else if(job$partial_info == 1) {
    rates <- calc_DE_discrepancy(data$simulation$abundances,
                                 data$simulation$observed_counts2,
                                 data$simulation$groups,
                                 method = method,
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
