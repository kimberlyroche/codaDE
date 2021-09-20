source("path_fix.R")

library(codaDE)
library(tidyverse)
library(RColorBrewer)
library(optparse)

DE_methods <- c("ALDEx2", "DESeq2", "scran")
use_baseline <- "self"

fit_predictive_model(DE_methods = DE_methods,
                     use_baseline = use_baseline,
                     output_weights = TRUE,
                     exclude_partials = TRUE,
                     exclude_independent = FALSE,
                     fit_full = TRUE)

for(use_result_type in c("TPR", "FPR")) {
  # Evaluate TPR model
  fit_obj <- readRDS(file.path("output",
                               "predictive_fits",
                               paste0(use_baseline, "_", use_result_type, ".rds")))

  test_data <- cbind(fit_obj$test_features, response = fit_obj$test_response)
  
  prediction <- predict(fit_obj$result, newdata = test_data) # predict.all = TRUE

  cat(paste0("R^2 (",use_result_type,"): ",
             round(cor(fit_obj$test_response, prediction)^2, 3), "\n"))
}
