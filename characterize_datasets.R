source("path_fix.R")

library(tidyverse)
library(codaDE)
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

output_dir <- file.path("output", "datasets")
results <- NULL

for(i in 1:nrow(wishlist)) {
  cat(paste0("Processing job ", i, " / ", nrow(wishlist), "\n"))

  # Parse data set
  job <- wishlist[i,]
  data <- readRDS(file.path(output_dir, paste0(job$uuid, ".rds")))

  abundances <- data$simulation$abundances
  if(job$partial_info == 1) {
    counts <- data$simulation$observed_counts2
  } else {
    counts <- data$simulation$observed_counts1
  }
  
  med_abs <- mean(rowSums(abundances))
  med_rel <- mean(rowSums(counts))

  n <- nrow(counts)/2
  counts_A <- counts[1:n,]
  counts_B <- counts[(n+1):(n*2),]
  
  results_row <- cbind(job, med_abs, med_rel, characterize_dataset(counts_A, counts_B))
  
  if(is.null(results)) {
    results <- results_row
  } else {
    results <- rbind(results, results_row)
  }
}

output_fn <- file.path("temp", paste0(output, "_", start, "-", end, ".txt"))
write.table(results, file = output_fn)
