source("path_fix.R")

library(tidyverse)
library(codaDE)
library(cowplot)
library(RColorBrewer)

source("ggplot_fix.R")

# Create directories manually
dir.create("output", showWarnings = FALSE)
dir.create(file.path("output", "images"), showWarnings = FALSE)

datasets <- c("VieiraSilva", "Barlow", "Song", "Monaco", "Hagai", "Owens", "Klein", "Yu")
thresholds <- c(1, 1, 1, 2, 1, 1, 1, 1)
# thresholds <- rep(2, length(datasets))
names(thresholds) <- datasets
model_dir <- "self_nopartial"

# ------------------------------------------------------------------------------
#   Classification accuracy
# ------------------------------------------------------------------------------

# Pull results files matching pattern

result_files <- list.files(path = file.path("output",
                                            "predictive_fits",
                                            model_dir,
                                            "classification",
                                            "validation_results",
                                            "no_norm"),
                           pattern = "results_(.*?)_threshold(\\d+)\\.tsv",
                           full.names = TRUE)
results <- NULL
for(file in result_files) {
  temp <- read.table(file, sep = "\t", header = TRUE)
  results <- rbind(results,
                   temp)
}

# Filter for the correct thresholds
keep_idx <- sapply(1:nrow(results), function(i) {
  results$threshold[i] == thresholds[[results$dataset[i]]]
})
results <- results[keep_idx,]

for(this_DE_method in unique(results$DE_method) ) {
  results_wrangled <- results %>%
    filter(DE_method == this_DE_method) %>%
    select(dataset, result_type, score_type, point) %>%
    pivot_wider(names_from = "result_type", values_from = "point") %>%
    mutate(agree = (true == predicted))
  tpr_agree_vec <- results_wrangled %>%
    filter(score_type == "TPR") %>%
    pull(agree)
  misses_tpr <- results_wrangled %>%
    filter(score_type == "TPR" & !agree) %>%
    pull(dataset)
  fpr_agree_vec <- results_wrangled %>%
    filter(score_type == "FPR") %>%
    pull(agree)
  misses_fpr <- results_wrangled %>%
    filter(score_type == "FPR" & !agree) %>%
    pull(dataset)
  cat(paste0(this_DE_method, " accuracy (TPR): ", sum(tpr_agree_vec) / length(tpr_agree_vec), " -- "))
  cat(paste0("missed on: ", paste0(misses_tpr, collapse = " "), "\n"))
  cat(paste0(this_DE_method, " accuracy (FPR): ", sum(fpr_agree_vec) / length(fpr_agree_vec), " -- "))
  cat(paste0("missed on: ", paste0(misses_fpr, collapse = " "), "\n"))
}

# ------------------------------------------------------------------------------
#   Regression accuracy
# ------------------------------------------------------------------------------

palette <- list(ALDEx2 = "#46A06B",
                DESeq2 = "#FF5733",
                # MAST = "#EF82BB",
                # NBGLM = "#7E54DE",
                scran = "#E3C012")

plot_labels <- list(FPR = "specificity (1 - FPR)", TPR = "sensitivity (TPR)")

# Pull results files matching pattern

result_files <- list.files(path = file.path("output",
                                            "predictive_fits",
                                            model_dir,
                                            "regression",
                                            "validation_results",
                                            "no_norm"),
                           pattern = "results_(.*?)_threshold(\\d+)\\.tsv",
                           full.names = TRUE)
results <- NULL
for(file in result_files) {
  temp <- read.table(file, sep = "\t", header = TRUE)
  results <- rbind(results,
                   temp)
}

# This loop also calculates R^2 for the real data x predictions
for(use_result_type in c("TPR", "FPR")) {
  plots <- list()
  legend <- NULL
  rsquared_df <- NULL
  for(this_dataset in datasets) {
    # Wrangle the poorly organized data
    plot_df <- results %>%
      filter(dataset == this_dataset) %>%
      filter(threshold == thresholds[[this_dataset]]) %>%
      filter(result_type == "predicted") %>%
      filter(score_type == use_result_type) %>%
      select(!c(threshold, score_type))
    plot_df$true <- NA
    for(i in 1:nrow(plot_df)) {
      plot_df$true[i] <- results %>%
        filter(dataset == this_dataset) %>%
        filter(threshold == thresholds[[this_dataset]]) %>%
        filter(result_type == "true") %>%
        filter(score_type == use_result_type) %>%
        filter(DE_method == plot_df$DE_method[i]) %>%
        pull(point)
    }
    rsquared_df <- rbind(rsquared_df,
                         plot_df %>% select(dataset, DE_method, point, true))

    pl <- ggplot(plot_df, aes(x = true, y = point)) +
      geom_segment(data = data.frame(x = 0, xend = 1, y = 0, yend = 1),
                   mapping = aes(x = x, xend = xend, y = y, yend = yend)) +
      geom_linerange(data = plot_df,
                     aes(ymin = lower90, ymax = upper90, color = factor(DE_method)),
                     size = 0.75) +
      geom_linerange(data = plot_df,
                     aes(ymin = lower50, ymax = upper50, color = factor(DE_method)),
                     size = 2) +
      geom_point(data = plot_df,
                 aes(x = true, y = point, color = factor(DE_method)),
                 size = 4) +
      scale_color_manual(values = palette) +
      theme_bw() +
      ylim(c(0,1)) +
      labs(x = paste0("observed ", plot_labels[[use_result_type]], "\n", this_dataset, " et al."),
           y = paste0("predicted ", plot_labels[[use_result_type]]),
           color = "Data type") +
      theme(axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)))
    legend <- get_legend(pl + theme(legend.position = "bottom"))
    pl <- pl +
      theme(legend.position = "none")
    
    plots[[this_dataset]] <- pl
  }
  
  # Report R^2
  cat(paste0("Overall ", use_result_type, " R^2: ", round(cor(rsquared_df$point, rsquared_df$true)**2, 3), "\n"))
  for(DE_method in sort(unique(rsquared_df$DE_method))) {
    idx <- which(rsquared_df$DE_method == DE_method)
    cat(paste0("for ", DE_method, ": ", round(cor(rsquared_df$point[idx], rsquared_df$true[idx])**2, 3), "\n"))
  }

  # Plot
  prow1 <- plot_grid(plotlist = plots[1:4], ncol = 4, labels = c("a", "b", "c", "d"),
                     #hjust = -1.35,
                     vjust = -0.25, label_size = 24)
  prow2 <- plot_grid(plotlist = plots[5:8], ncol = 4, labels = c("e", "f", "g", "h"),
                     #hjust = -1.35,
                     vjust = -0.2, label_size = 24)
  pl <- plot_grid(prow1, prow2, legend, ncol = 1, rel_heights = c(1, 1, .1), scale = 0.9)
  pl <- pl +
    theme(plot.margin = ggplot2::margin(t = 25, r = 0, b = 0, l = 0, "pt"))
  show(pl)
  # ggsave(file.path("output",
  #                  "images",
  #                  paste0("validations_",
  #                         use_result_type,
  #                         ".png")),
  #        plot = pl,
  #        dpi = 100,
  #        units = "in",
  #        height = 7.35,
  #        width = 13)
}
