# This code runs and evaluates error in differential expression calls for a single
# simulation instance with the desired characteristics.

library(codaDE)
library(optparse)

option_list = list(
  make_option(c("--p"), type = "numeric", default = NULL, 
              help = "number of genes", metavar = "numeric"),
  make_option(c("--n"), type = "numeric", default = NULL, 
              help = "number of samples per condition", metavar = "numeric"),
  make_option(c("--k"), type = "numeric", default = NULL, 
              help = "number of cell types (if simulating single-cell)", metavar = "numeric"),
  make_option(c("--run_label"), type = "character", default = NULL, 
              help = "label associated with this simulation (to appear in output file)", metavar = "character"),
  make_option(c("--NB_for_DE"), type = "logical", default = TRUE,
              help = "flag associated with method to use for differential expression calling: TRUE = NB + LRT, FALSE = log-LM + permutations", metavar = "numeric"),
  make_option(c("--use_ALR"), type = "logical", default = FALSE, 
              help = "flag indicating whether or not to perform differential abundance testing on logratios", metavar = "logical"),
  make_option(c("--filter_abundance"), type = "numeric", default = 0, 
              help = "minimum average abundance of a gene on which to test differential abundance", metavar = "numeric"),
  make_option(c("--rarefy"), type = "logical", default = FALSE, 
              help = "flag indicating whether or not to rarefy samples", metavar = "logical")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

sweep_simulations(p = opt$p,
                  n = opt$n,
                  run_label = opt$run_label,
                  k = opt$k,
                  de_sweep = c((1/10), (1/4), (2/3)),
                  corr_sweep = c(0, 0.5, 0.9),
                  output_file = paste0("results_",opt$run_label,".tsv"),
                  alpha = 0.05,
                  use_ALR = opt$use_ALR,
                  filter_abundance = opt$filter_abundance,
                  call_DA_by_NB = opt$NB_for_DE,
                  rarefy = opt$rarefy)
