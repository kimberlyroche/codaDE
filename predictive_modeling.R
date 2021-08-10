source("path_fix.R")

library(codaDE)
library(tidyverse)
library(RColorBrewer)

use_methods <- c("all")
# use_methods <- c("ALDEx2", "DESeq2", "MAST", "scran")

plot_labels <- list(fpr = "specificity (1 - FPR)", tpr = "sensitivity (TPR)")
save_tag <- paste0("_", as.numeric(Sys.time()))

for(DE_method in use_methods) {
#  for(use_baseline in c("self", "oracle")) {
  for(use_baseline in c("self")) {
    for(use_result_type in c("TPR", "FPR")) {
      save_fn <- file.path("output",
                           "predictive_fits",
                           DE_method,
                           paste0(DE_method, "_", use_result_type, "_", use_baseline, save_tag, ".rds"))
      if(!file.exists(save_fn)) {
        cat(paste0("Fitting model for ", DE_method, " on ", use_result_type, " (", use_baseline, ")\n"))
        fit_predictive_model(DE_method = DE_method,
                             plot_weights = TRUE,
                             exclude_partials = TRUE,
                             exclude_independent = FALSE,
                             train_percent = 0.5,
                             save_tag = save_tag)
      }

      fit_obj <- readRDS(save_fn)
      
      rf_test_data <- cbind(fit_obj$test_features, fit_obj$test_response)
      prediction <- predict(fit_obj$result, newdata = rf_test_data)
      
      plot_df <- data.frame(true = fit_obj$test_response,
                            predicted = prediction,
                            p = fit_obj$test_features$P)
      
      # Alternative, 2D density plot
      # pl <- ggplot(plot_df, aes(x = true, y = predicted)) +
      #   geom_density_2d(color = "black") +
      #   xlim(c(0,1)) +
      #   ylim(c(0,1)) +
      #   labs(x = paste0("observed ", plot_labels[[use_result_type]]),
      #        y = paste0("predicted ", plot_labels[[use_result_type]]),
      #        fill = "feature no.")
      
      # pl <- ggplot(plot_df, aes(x = true, y = predicted, fill = factor(p))) +
      #   geom_point(shape = 21, size = 3) +
      #   scale_fill_brewer(palette = "RdYlGn") +
      #   xlim(c(0,1)) +
      #   ylim(c(0,1)) +
      #   labs(x = paste0("observed ", plot_labels[[use_result_type]]),
      #        y = paste0("predicted ", plot_labels[[use_result_type]]),
      #        fill = "Feature number")
      # show(pl)
      # ggsave(file.path("output",
      #                  "predictive_fits",
      #                  DE_method,
      #                  paste0("predictions_",
      #                         DE_method,
      #                         "_",
      #                         use_result_type,
      #                         "_",
      #                         use_baseline,
      #                         ".png")),
      #        plot = pl,
      #        dpi = 100,
      #        units = "in",
      #        height = 4,
      #        width = 5.5)
      
      cat(paste0("R^2 (",toupper(use_result_type),"): ",
                 round(cor(plot_df$true, plot_df$predicted)^2, 3), "\n"))
    }
  }
}
