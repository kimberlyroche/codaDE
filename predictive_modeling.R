source("path_fix.R")

library(codaDE)
library(tidyverse)
library(RColorBrewer)

DE_method <- "all"
use_baseline <- "self"

plot_labels <- list(fpr = "specificity (1 - FPR)", tpr = "sensitivity (TPR)")

# Append as.numeric(Sys.time()) if (for example) training many models on subsets
# of the data, in order to make model fit file names unique
save_slug <- file.path("output",
                       "predictive_fits",
                       DE_method,
                       paste0(DE_method, "_", use_baseline))
fit_predictive_model(DE_method = DE_method,
                     use_baseline = use_baseline,
                     plot_weights = FALSE,
                     exclude_partials = FALSE,
                     exclude_independent = FALSE)

use_result_type <- "TPR"

for(use_result_type in c("TPR", "FPR")) {
  # Evaluate TPR model
  fit_obj <- readRDS(paste0(save_slug, "_", use_result_type, ".rds"))
  
  str(fit_obj, max.level = 1)
      
  rf_test_data <- cbind(fit_obj$test_features, response = fit_obj$test_response)
  head(rf_test_data)
  
  prediction <- predict(fit_obj$result, newdata = rf_test_data, predict.all = TRUE)

  plot_df <- data.frame(true = fit_obj$test_response,
                        predicted = prediction$aggregate,
                        p = fit_obj$test_features$P)
  plot_df$lower <- apply(prediction$individual, 1, function(x) quantile(x, probs = 0.05))
  plot_df$upper <- apply(prediction$individual, 1, function(x) quantile(x, probs = 0.95))

  cat(paste0("R^2 (",use_result_type,"): ",
             round(cor(plot_df$true, plot_df$predicted)^2, 3), "\n"))
  
  plot_df <- plot_df %>%
    arrange(true)
  ggplot(plot_df, aes(x = 1:nrow(plot_df), ymin = lower, y = true, ymax = upper)) +
    geom_ribbon(fill = "grey70") +
    geom_line() +
    labs(x = "sample index",
         y = paste0("prediction interval and true ", use_result_type))
}
