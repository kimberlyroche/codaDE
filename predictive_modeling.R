source("path_fix.R")

library(codaDE)
library(tidyverse)
library(caret)
library(RColorBrewer)
library(optparse)

DE_methods <- c("ALDEx2", "DESeq2", "scran")
use_baseline <- "self"

fit_predictive_model(DE_methods = DE_methods,
                     use_baseline = use_baseline,
                     output_weights = TRUE,
                     exclude_partials = TRUE,
                     exclude_independent = FALSE,
                     do_classify = TRUE)
                     # abs_feature_list = c("FC_ABSOLUTE"))

for(use_result_type in c("TPR", "FPR")) {
  # Evaluate TPR model
  fit_obj <- readRDS(file.path("output",
                               "predictive_fits",
                               paste0(use_baseline, "_", use_result_type, ".rds")))

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
      acc_str <- paste0(prepend_str, "R^2 (",use_result_type,"): ",
                        round(cor(fit_obj$test_response, prediction)^2, 3))
    } else {
      # Classification
      acc_obj <- confusionMatrix(test_data$response[idx], prediction)
      acc_str <- paste0(prepend_str, "accuracy: ", round(acc_obj$overall[["Accuracy"]], 2), " (",
                        "sensitivity: ", round(acc_obj$byClass[["Sensitivity"]], 2), " / ",
                        "specificity: ", round(acc_obj$byClass[["Specificity"]], 2), ")")
    }

    cat(paste0(acc_str, "\n"))
  }
}
