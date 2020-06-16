# This code runs and evaluates error in differential expression calls for a single simulation instance with
# the desired characteristics.

library(codaDE)

args <- commandArgs(trailing = TRUE)
if(length(args) < 4) {
  stop("Arguments: {p} {fold change} {prop. DE features} {run label} {opt: mean size factor} {opt: size factor correlation} {opt: bimodal = T/F}")
}
p <- as.numeric(args[1])
if(p < 2) {
  stop("p must be at least 2!\n")
}
fc <- as.numeric(args[2])
if(fc < 0) {
  stop("Fold change must be greater than zero!\n")
}
ppf <- as.numeric(args[3])
if(ppf > 1 | ppf < 0) {
  stop("Proportion of perturbed features must be in range [0,1]!\n")
}
run_label <- args[4]
if(length(args) >= 5) {
  msf <- as.numeric(args[5])
  if(msf < 100) {
    stop("Mean size factor must be at least 100!\n")
  }
} else {
  msf <- 5000
}
if(length(args) >= 6) {
  sfc <- as.numeric(args[6])
  if(sfc > 1 | sfc < 0) {
    stop("Size factor correlation must be in range [0, 1]!\n")
  }
} else {
  sfc <- 0
}
if(length(args) >= 7) {
  bimodal <- as.logical(args[7])
} else {
  bimodal <- FALSE
}

run_evaluation_instance(p = p, fold_change = fc, proportion_perturbed_features = ppf, run_label = run_label, mean_size_factor = msf,
                        size_factor_correlation = sfc, bimodal = bimodal, output_file = "results.txt", alpha = 0.05)
