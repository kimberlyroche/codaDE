source("path_fix.R")

library(codaDE)
library(tidyverse)
library(RColorBrewer)

DE_method <- "all"
use_baseline <- "self"

plot_labels <- list(fpr = "specificity (1 - FPR)", tpr = "sensitivity (TPR)")

save_slug <- file.path("output",
                       "predictive_fits",
                       DE_method,
                       paste0(DE_method, "_", use_baseline, "_", as.numeric(Sys.time())))
fit_predictive_model(DE_method = DE_method,
                     use_baseline = use_baseline,
                     plot_weights = TRUE,
                     exclude_partials = TRUE,
                     exclude_independent = FALSE,
                     train_percent = 0.2,
                     save_slug = save_slug,
                     save_training_data = FALSE)

for(use_result_type in c("TPR", "FPR")) {
  # Evaluate TPR model
  fit_obj <- readRDS(paste0(save_slug, "_", use_result_type, ".rds"))
      
  rf_test_data <- cbind(fit_obj$test_features, fit_obj$test_response)
  prediction <- predict(fit_obj$result, newdata = rf_test_data)

  plot_df <- data.frame(true = fit_obj$test_response,
                        predicted = prediction,
                        p = fit_obj$test_features$P)
      
  cat(paste0("R^2 (",use_result_type,"): ",
             round(cor(plot_df$true, plot_df$predicted)^2, 3), "\n"))
}
