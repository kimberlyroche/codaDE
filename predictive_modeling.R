source("path_fix.R")

library(codaDE)
library(tidyverse)
library(glmnet)
library(caret)
library(RColorBrewer)
library(optparse)

option_list = list(
  make_option(c("--selfbaseline"),
              type = "logical",
              default = "FALSE",
              help = "use self calls as reference instead of oracle",
              metavar = "logical"),
  make_option(c("--usetotals"),
              type = "logical",
              default = "FALSE",
              help = "use information about change in total abundances from absolute counts",
              metavar = "logical"),
  make_option(c("--userenormcounts"),
              type = "logical",
              default = "FALSE",
              help = "use features generated from renormalized counts",
              metavar = "logical"),
  make_option(c("--cpm"),
              type = "logical",
              default = "FALSE",
              help = "convert relative abundances to counts per million",
              metavar = "logical")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

use_self_baseline <- opt$selfbaseline
use_totals <- opt$usetotals
use_renorm_counts <- opt$userenormcounts
use_cpm <- opt$cpm

DE_methods <- c("ALDEx2", "DESeq2", "scran")

model_dir <- ifelse(use_self_baseline, "self", "oracle")

for(DE_method in DE_methods) {
  fit_predictive_model(DE_methods = c(DE_method),
                       use_baseline = ifelse(use_self_baseline, "self", "oracle"),
                       use_totals = use_totals,
                       use_renorm_counts = use_renorm_counts,
                       output_weights = TRUE,
                       train_percent = 0.8,
                       use_cpm = use_cpm)
}

# Evaluate performance
for(use_result_type in c("TPR", "FPR")) {
  for(DE_method in DE_methods) {
    save_fn <- file.path("output",
                         "predictive_fits",
                         model_dir,
                         "regression",
                         paste0(use_result_type,
                                "_",
                                DE_method,
                                ".rds"))
    if(file.exists(save_fn)) {
      fit_obj <- readRDS(save_fn)
      test_data <- cbind(fit_obj$test_features, response = fit_obj$test_response)
      prediction <- predict(fit_obj$result, newdata = fit_obj$test_features) # predict.all = TRUE
      
      acc_str <- paste0(DE_method, " R^2 (", use_result_type, "): ",
                        round(cor(fit_obj$test_response, prediction)^2, 3))
      
      cat(paste0(acc_str, "\n"))
    }
  }
}
