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

output_fn <- file.path("temp", paste0(output, "_", start, "-", end, ".txt"))

# Initialize output
output_file <- file(output_fn)
writeLines(paste0(c("ID", "uuid", "measured_by", "percent_diff_realiz"),
                  collapse = "\t"),
           output_file)
close(output_file)

# ------------------------------------------------------------------------------
#  Pull jobs to run, then run 'em!
# ------------------------------------------------------------------------------

wishlist <- read.table(input)
start <- min(start, nrow(wishlist))
end <- min(end, nrow(wishlist))
wishlist <- wishlist[start:end,]

output_dir <- file.path("output", "datasets")
results <- NULL

counter <- 1
for(i in 1:nrow(wishlist)) {
  cat(paste0("Processing job ", i, " / ", nrow(wishlist), "\n"))

  # Parse data set
  job <- wishlist[i,]
  data <- readRDS(file.path(output_dir, paste0(job$uuid, ".rds")))

  abundances <- data$simulation$abundances
  counts <- data$simulation$observed_counts1
  
  results_row <- cbind(job, "fc_1.5", calc_threshold_DA(counts, fc_lower = 1/1.5, fc_upper = 1.5))
  
  write_delim(results_row, output_fn, delim = "\t", append = TRUE)
  counter <- counter + 1
}
