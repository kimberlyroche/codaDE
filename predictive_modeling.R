source("path_fix.R")

library(codaDE)
library(tidyverse)
library(optparse)

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

plot_density <- FALSE
plot_labels <- list(fpr = "specificity (1 - FPR)", tpr = "sensitivity (TPR)")

# Fit predictive model (w/ predictions!)
for(DE_method in use_methods) {
  fit_obj <- fit_predictive_model(model = model, do_predict = TRUE,
                                  DE_method = DE_method, plot_weights = TRUE)

  for(use_result_type in c("tpr", "fpr")) {
    test_response <- fit_obj$predictions[[use_result_type]]$true
    prediction <- fit_obj$predictions[[use_result_type]]$predicted
    p_labels <- fit_obj$predictions[[use_result_type]]$p_labels
    
    plot_df <- data.frame(true = test_response,
                          predicted = prediction,
                          p = factor(p_labels))
    levels(plot_df$p) <- c("low", "med", "high")
    if(plot_density) {
      pl <- ggplot(plot_df, aes(x = true, y = predicted)) +
        geom_density_2d(color = "black") +
        xlim(c(0,1)) +
        ylim(c(0,1)) +
        labs(x = paste0("observed ", plot_labels[[use_result_type]]),
             y = paste0("predicted ", plot_labels[[use_result_type]]),
             fill = "feature no.")
    } else {
      pl <- ggplot(plot_df, aes(x = true, y = predicted, fill = p)) +
        geom_point(shape = 21, size = 3) +
        xlim(c(0,1)) +
        ylim(c(0,1)) +
        labs(x = paste0("observed ", plot_labels[[use_result_type]]),
             y = paste0("predicted ", plot_labels[[use_result_type]]),
             fill = "feature no.")
    }
    
    show(pl)
    ggsave(file.path("output",
                     "images",
                     paste0(model, "_results"),
                     paste0(model,
                            "_predictions_",
                            DE_method,
                            "_",
                            use_result_type,
                            ".png")),
           plot = pl,
           dpi = 100,
           units = "in",
           height = 4,
           width = 5)
    
    cat(paste0("R^2 (",toupper(use_result_type),"): ",
               round(cor(test_response, prediction)^2, 3), "\n"))
  }
}

