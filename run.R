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
  make_option(c("--method"), type = "character", default = NULL, 
              help = "differential abundance calling method", metavar = "character"),
  make_option(c("--use_ALR"), type = "logical", default = FALSE, 
              help = "flag indicating whether or not to perform differential abundance testing on logratios", metavar = "logical"),
  make_option(c("--filter_abundance"), type = "numeric", default = 0, 
              help = "minimum average abundance of a gene on which to test differential abundance", metavar = "numeric")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

de_sweep <- c(0.2, 0.4, 0.6, 0.8)
if(!is.null(opt$prop_da)) {
  de_sweep <- c(opt$prop_da)
}

corr_sweep <- c(0, 0.5, 0.9)
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
                  methods = c(opt$method))
