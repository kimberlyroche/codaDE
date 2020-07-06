# This code runs and evaluates error in differential expression calls for a single simulation instance with
# the desired characteristics.

library(codaDE)

args <- commandArgs(trailing = TRUE)
if(length(args) < 2) {
  stop("Missing arguments: {number of genes} {number of samples per treatment} {run label}")
}
p <- as.numeric(args[1])
n <- as.numeric(args[2])
run_label <- args[3]
# optional arguments
use_ALR <- FALSE
if(length(args) > 3) {
  use_ALR <- as.logical(args[4])
}
filter_abundance <- 0
if(length(args) > 4) {
  filter_abundance <- as.numeric(args[5])
}
rarefy <- FALSE
if(length(args) > 5) {
  rarefy <- as.logical(args[6])
}

sweep_RNAseq(p,
             n,
             run_label = run_label,
             de_sweep = c((1/10), (1/4), (2/3)),
             corr_sweep = c(0, 0.9),
             output_file = paste0("results_",run_label,".tsv"),
             use_ALR = use_ALR,
             filter_abundance = filter_abundance,
             rarefy = rarefy)
