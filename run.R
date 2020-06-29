# This code runs and evaluates error in differential expression calls for a single simulation instance with
# the desired characteristics.

library(codaDE)

args <- commandArgs(trailing = TRUE)
if(length(args) < 2) {
  stop("Missing arguments: {number of genes} {number of samples per treatment}")
}
p <- as.numeric(args[1])
n <- as.numeric(args[2])
rarefy <- FALSE
if(length(args) > 2) {
  rarefy <- as.logical(args[3])
}
use_ALR <- FALSE
if(length(args) > 3) {
  use_ALR <- as.logical(args[4])
}

run_label <- "RNAseq_like_CODA"

sweep_RNAseq(p, n, run_label = run_label, de_sweep = seq(from = 0.1, to = 0.9, by = 0.1), corr_sweep = seq(from = 0.1, to = 0.9, by = 0.1),
             output_file = "results.txt", alpha = 0.05, use_ALR = use_ALR, rarefy = rarefy)
