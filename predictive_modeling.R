source("path_fix.R")

library(codaDE)
library(tidyverse)
library(glmnet)
library(caret)
library(RColorBrewer)
library(optparse)

option_list = list(
  make_option(c("--classify"),
              type = "logical",
              default = "FALSE",
              help = "classification (vs. regression) flag",
              metavar = "logical"),
  make_option(c("--alpha"),
              type = "numeric",
              default = "0.95",
              help = "classification threshold",
              metavar = "numeric"),
  make_option(c("--permodel"),
              type = "logical",
              default = "TRUE",
              help = "use per-model predictive fits",
              metavar = "logical"),
  make_option(c("--selfbaseline"),
              type = "logical",
              default = "TRUE",
              help = "use self calls as reference instead of oracle",
              metavar = "logical"),
  make_option(c("--partials"),
              type = "logical",
              default = "FALSE",
              help = "flag indicating whether or not to use simulations with partially informative total abundances",
              metavar = "logical")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

do_classify <- opt$classify
alpha <- opt$alpha
per_model <- opt$permodel
use_self_baseline <- opt$selfbaseline
use_partials <- opt$partials

model_type <- "RF"

DE_methods <- c("ALDEx2", "DESeq2", "scran")

model_dir <- paste0(ifelse(use_self_baseline, "self", "oracle"),
                    "_",
                    ifelse(use_partials, "partial", "nopartial"))

for(use_result_type in c("TPR", "FPR")) {
  if(per_model) {
    for(DE_method in DE_methods) {
      save_fn <- file.path("output",
                           "predictive_fits",
                           model_dir,
                           ifelse(do_classify, "classification", "regression"),
                                paste0(use_result_type,
                                          "_",
                                          DE_method,
                                          ".rds"))
      if(!file.exists(save_fn)) {
        fit_predictive_model(model_type = model_type,
                             DE_methods = DE_method,
                             use_baseline = ifelse(use_self_baseline, "self", "oracle"),
                             output_weights = TRUE,
                             exclude_partials = !use_partials,
                             exclude_independent = FALSE,
                             do_classify = do_classify,
                             alpha = alpha)
      }
      fit_obj <- readRDS(save_fn)

      test_data <- cbind(fit_obj$test_features, response = fit_obj$test_response)
      
      if(model_type == "LR") {
        newx <- model.matrix(response ~ ., test_data)
        prediction <- factor(predict(fit_obj$result, newx = newx, s = "lambda.1se", type = "class"))
      } else {
        prediction <- predict(fit_obj$result, newx = test_data) # predict.all = TRUE
      }

      if(any(!(test_data$response %in% c(0,1)))) {
        # Regression
        acc_str <- paste0(DE_method, " R^2 (", use_result_type, "): ",
                          round(cor(fit_obj$test_response, prediction)^2, 3))
      } else {
        # Classification
        acc_obj <- confusionMatrix(test_data$response, prediction)
        acc_str <- paste0(DE_method, " accuracy (", use_result_type, "): ", round(acc_obj$overall[["Accuracy"]], 2), " (",
                          "sensitivity: ", round(acc_obj$byClass[["Sensitivity"]], 2), " / ",
                          "specificity: ", round(acc_obj$byClass[["Specificity"]], 2), ")")
      }

      cat(paste0(acc_str, "\n"))
    }
  } else {
    save_fn <- file.path("output",
                         "predictive_fits",
                         model_dir,
                         ifelse(do_classify, "classification", "regression"),
                         paste0(use_result_type, ".rds"))
    if(!file.exists(save_fn)) {
      fit_predictive_model(model_type = model_type,
                           DE_methods = DE_methods,
                           use_baseline = ifelse(use_self_baseline, "self", "oracle"),
                           output_weights = TRUE,
                           exclude_partials = !use_partials,
                           exclude_independent = FALSE,
                           do_classify = do_classify,
                           alpha = alpha)
    }
    fit_obj <- readRDS(save_fn)

    test_data <- cbind(fit_obj$test_features, response = fit_obj$test_response)

    for(p in c(-1, 100, 1000, 5000)) {
      if(p > 0) {
        idx <- which(test_data$P == p)
      } else {
        idx <- 1:nrow(test_data)
      }
      prediction <- predict(fit_obj$result, newdata = test_data[idx,]) # predict.all = TRUE

      prepend_str <- "Overall "
      if(p > 0) {
        prepend_str <- paste0("P=", p, " ")
      }

      if(any(!(test_data$response %in% c(0,1)))) {
        # Regression
        acc_str <- paste0(prepend_str, "R^2 (", use_result_type, "): ",
                          round(cor(fit_obj$test_response, prediction)^2, 3))
      } else {
        # Classification
        acc_obj <- confusionMatrix(test_data$response[idx], prediction)
        acc_str <- paste0(prepend_str, "accuracy (", use_result_type, "): ", round(acc_obj$overall[["Accuracy"]], 2), " (",
                          "sensitivity: ", round(acc_obj$byClass[["Sensitivity"]], 2), " / ",
                          "specificity: ", round(acc_obj$byClass[["Specificity"]], 2), ")")
      }

      cat(paste0(acc_str, "\n"))
    }
  }
}
