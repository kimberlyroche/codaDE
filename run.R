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
  make_option(c("--prop_da"), type = "numeric", default = NULL, 
              help = "proportion of differential expression", metavar = "numeric"),
  make_option(c("--sf_corr"), type = "numeric", default = NULL, 
              help = "size factor correlation", metavar = "numeric"),
  make_option(c("--NB_for_DE"), type = "logical", default = TRUE,
              help = "flag associated with method to use for differential expression calling: TRUE = NB + LRT, FALSE = log-LM + permutations", metavar = "numeric"),
  make_option(c("--use_ALR"), type = "logical", default = FALSE, 
              help = "flag indicating whether or not to perform differential abundance testing on logratios", metavar = "logical"),
  make_option(c("--filter_abundance"), type = "numeric", default = 0, 
              help = "minimum average abundance of a gene on which to test differential abundance", metavar = "numeric"),
  make_option(c("--rarefy"), type = "logical", default = FALSE, 
              help = "flag indicating whether or not to rarefy samples", metavar = "logical"),
  make_option(c("--existing"), type = "logical", default = FALSE, 
              help = "use existing simulation", metavar = "logical"),
  make_option(c("--save_slot"), type = "numeric", default = NULL, 
              help = "index in results to save evaluation run into", metavar = "numeric")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

de_sweep <- seq(from = 0.1, to = 0.9, by = 0.1)
de_sweep <- c(0.3, 0.7, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 0.9)
if(!is.null(opt$prop_da)) {
  de_sweep <- c(opt$prop_da)
}

corr_sweep <- seq(from = 0.1, to = 0.9, by = 0.1)
corr_sweep <- c(0.3, 0.7, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 0.9)
if(!is.null(opt$sf_corr)) {
  corr_sweep <- c(opt$sf_corr)
}

sweep_simulations(p = opt$p,
                  n = opt$n,
                  k = opt$k,
                  de_sweep = de_sweep,
                  corr_sweep = corr_sweep,
                  use_ALR = opt$use_ALR,
                  filter_abundance = opt$filter_abundance,
                  call_DA_by_NB = opt$NB_for_DE,
                  rarefy = opt$rarefy,
                  use_existing_simulations = opt$existing,
                  save_slot = opt$save_slot)
