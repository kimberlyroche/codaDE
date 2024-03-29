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
                    "type",
                    "partial_info",
                    "method",
                    "tpr",
                    "fpr",
                    "fdr",
                    "beta"), collapse = "\t"),
           output_file)
close(output_file)

wishlist <- read.table(input)
start <- min(start, nrow(wishlist))
end <- min(end, nrow(wishlist))
wishlist <- wishlist[start:end,]

counter <- 1
for(i in 1:nrow(wishlist)) {
  job <- wishlist[i,]
  all_calls <- list(calls = list(pval = job$calls, betas = job$betas),
                    oracle_calls = list(pval = job$baseline_calls, betas = job$baseline_betas))
  rates <- calc_DA_discrepancy(all_calls$calls, all_calls$oracle_calls, alpha = job$fdr, beta = job$beta)
  results_row <- data.frame(ID = counter,
                            uuid = job$uuid,
                            baseline = job$baseline,
                            type = job$type,
                            partial_info = job$partial_info,
                            method = job$method,
                            tpr = rates$TPR,
                            fpr = rates$FPR,
                            fdr = job$fdr,
                            beta = job$beta)
  write_delim(results_row, output_fn, delim = "\t", append = TRUE)
  counter <- counter + 1
}
cat(paste0(counter, " rows updated\n"))
