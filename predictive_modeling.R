source("path_fix.R")

library(codaDE)
library(tidyverse)
library(RColorBrewer)
library(optparse)

DE_method <- "all"
use_baseline <- "self"
model_type <- "RF"

plot_labels <- list(fpr = "specificity (1 - FPR)", tpr = "sensitivity (TPR)")

# Append as.numeric(Sys.time()) to `save_slug` if (for example) training many 
# models on subsets of the data, in order to make model fit file names unique
fit_predictive_model(model_type = model_type,
                     DE_method = DE_method,
                     use_baseline = use_baseline,
                     plot_weights = FALSE,
                     exclude_partials = FALSE,
                     exclude_independent = FALSE)

save_slug <- file.path("output",
                       "predictive_fits",
                       DE_method,
                       paste0(DE_method, "_", use_baseline))

for(use_result_type in c("TPR", "FPR")) {
  # Evaluate TPR model
  fit_obj <- readRDS(paste0(save_slug, "_", model_type, "_", use_result_type, ".rds"))
      
  test_data <- cbind(fit_obj$test_features, response = fit_obj$test_response)
  
  prediction <- predict(fit_obj$result, newdata = test_data) # predict.all = TRUE

  plot_df <- data.frame(true = fit_obj$test_response,
                        # predicted = prediction$aggregate,
                        predicted = prediction,
                        p = fit_obj$test_features$P)
  # plot_df$p05 <- apply(prediction$individual, 1, function(x) quantile(x, probs = 0.05))
  # plot_df$p25 <- apply(prediction$individual, 1, function(x) quantile(x, probs = 0.25))
  # plot_df$p75 <- apply(prediction$individual, 1, function(x) quantile(x, probs = 0.75))
  # plot_df$p95 <- apply(prediction$individual, 1, function(x) quantile(x, probs = 0.95))
  
  cat(paste0("R^2 (",use_result_type,"): ",
             round(cor(plot_df$true, plot_df$predicted)^2, 3), "\n"))
  
  # plot_df <- plot_df %>%
  #   arrange(true)
  # pl <- ggplot() +
  #   geom_ribbon(data = plot_df,
  #               mapping = aes(x = 1:nrow(plot_df), ymin = p05, ymax = p95),
  #               fill = "grey70") +
  #   geom_ribbon(data = plot_df,
  #               mapping = aes(x = 1:nrow(plot_df), ymin = p25, ymax = p75),
  #               fill = "grey50") +
  #   geom_point(data = plot_df,
  #              mapping = aes(x = 1:nrow(plot_df), y = true),
  #              size = 2) +
  #   labs(x = "sample index",
  #        y = paste0("prediction interval and true ", use_result_type)) +
  #   theme_bw()
  # show(pl)
  # ggsave(file.path(save_dir,
  #                  paste0("intervals_",
  #                         DE_method,
  #                         "_",
  #                         use_baseline,
  #                         "_",
  #                         use_result_type,
  #                         ".png")),
  #        plot = pl,
  #        dpi = 100,
  #        units = "in",
  #        height = 4,
  #        width = 6)
}





