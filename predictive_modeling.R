source("path_fix.R")

library(codaDE)
library(tidyverse)
library(optparse)
library(RColorBrewer)

option_list = list(
  make_option(c("--model"), type = "character", default = "EN",
              help = "predictive model to use: EN, GP, linear, RF", metavar = "character"),
  make_option(c("--together"), type = "logical", default = "FALSE",
              help = "combine results of all methods", metavar = "logical")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);
model <- opt$model
if(!(model %in% c("EN", "GP", "linear", "RF"))) {
  stop("Invalid model specified!\n")
}
together <- opt$together

model <- "RF"
if(together) {
  use_methods <- c("all")
} else {
  use_methods <- c("ALDEx2", "DESeq2", "MAST", "NBGLM", "scran")
}

ramped_palette <- colorRampPalette(c("#7479c4", "#e63030"))(3)

plot_density <- FALSE
plot_labels <- list(fpr = "specificity (1 - FPR)", tpr = "sensitivity (TPR)")

for(DE_method in use_methods) {
  fit_obj <- fit_predictive_model(model = model, do_predict = TRUE,
                                  DE_method = DE_method, plot_weights = TRUE,
                                  exclude_partials = TRUE, exclude_indepedent = FALSE)

  for(use_baseline in c("self", "threshold")) {
    for(use_result_type in c("tpr", "fpr")) {
      # Save fitted model -- already done by fit_predictive_model()
      # base_path <- c("output",
      #                "images",
      #                paste0(model, "_results"),
      #                ifelse(together, "all", DE_method))
      # for(i in 1:length(base_path)) {
      #   suppressWarnings(dir.create(do.call(file.path, as.list(base_path[1:i]))))
      # }
      # model_fn <- file.path(do.call(file.path, as.list(base_path)),
      #                       paste0(model,
      #                              "_",
      #                              ifelse(together, "all", DE_method),
      #                              "_",
      #                              use_result_type,
      #                              "_",
      #                              use_baseline,
      #                              ".rds"))
      # saveRDS(fit_obj$fitted_model[[use_baseline]][[use_result_type]], model_fn)
      
      test_response <- fit_obj$predictions[[use_baseline]][[use_result_type]]$true
      prediction <- fit_obj$predictions[[use_baseline]][[use_result_type]]$predicted
      p_labels <- fit_obj$predictions[[use_baseline]][[use_result_type]]$p_labels
      
      plot_df <- data.frame(true = test_response,
                            predicted = prediction,
                            p = p_labels)
      if(plot_density) {
        pl <- ggplot(plot_df, aes(x = true, y = predicted)) +
          geom_density_2d(color = "black") +
          xlim(c(0,1)) +
          ylim(c(0,1)) +
          labs(x = paste0("observed ", plot_labels[[use_result_type]]),
               y = paste0("predicted ", plot_labels[[use_result_type]]),
               fill = "feature no.")
      } else {
        pl <- ggplot(plot_df, aes(x = true, y = predicted, fill = factor(p))) +
          geom_point(shape = 21, size = 3) +
          scale_fill_brewer(palette = "RdYlGn") +
          xlim(c(0,1)) +
          ylim(c(0,1)) +
          labs(x = paste0("observed ", plot_labels[[use_result_type]]),
               y = paste0("predicted ", plot_labels[[use_result_type]]),
               fill = "Feature number")
      }
      
      show(pl)
      ggsave(file.path("output",
                       "images",
                       paste0(model, "_results"),
                       DE_method,
                       paste0(model,
                              "_predictions_",
                              DE_method,
                              "_",
                              use_result_type,
                              "_",
                              use_baseline,
                              ".png")),
             plot = pl,
             dpi = 100,
             units = "in",
             height = 4,
             width = 5.5)
      
      cat(paste0("R^2 (",toupper(use_result_type),"): ",
                 round(cor(test_response, prediction)^2, 3), "\n"))
    }
  }
}
