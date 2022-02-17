source("path_fix.R")

library(tidyverse)
library(codaDE)
library(cowplot)
library(RColorBrewer)

# Create directories manually
dir.create("output", showWarnings = FALSE)
dir.create(file.path("output", "images"), showWarnings = FALSE)

datasets <- c("VieiraSilva", "Muraro", "Hagai", "Hashimshony", "Gruen", "Kimmerling",
              "Song", "Barlow", "Monaco", "Yu", "Klein", "Owens")
thresholds <- rep(1, length(datasets))

names(thresholds) <- datasets
model_dir <- "oracle"
submodel_dir <- "regression_cpm"

# palette <- list(ALDEx2 = "#46A06B",
#                 DESeq2 = "#FF5733",
#                 scran = "#E3C012")

plot_labels <- list(FPR = "specificity (1 - FPR)", TPR = "sensitivity (TPR)")

result_files <- list.files(path = file.path("output",
                                            "predictive_fits",
                                            model_dir,
                                            submodel_dir,
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

# Filter out scran on Monaco et al. where the low sample number (4) causes
# the computeSumFactors() scaling factor estimate to throw an error. In general,
# scran seems to work best with an absolute MINIMUM of 5 samples.
results <- results %>%
  filter(!(dataset == "Monaco" & DE_method == "scran"))

plotlist <- list()
for(i in 1:length(datasets)) {
  this_dataset <- datasets[i]
  
  # ----------------------------------------------------------------------------
  #   Plot 1: Parse absolute counts of true positives, false positives, etc.
  # ----------------------------------------------------------------------------
  
  bar_components <- NULL
  for(method in c("ALDEx2", "DESeq2", "scran")) {
    if(this_dataset == "Monaco" & method == "scran") next;
    calls <- readRDS(file.path("output",
                               "real_data_calls",
                               "no_norm",
                               paste0("calls_oracle_",method,"_",this_dataset,"_threshold",thresholds[i],".rds")))
    TP <- sum(calls$rates$TP_calls)
    TN <- sum(calls$rates$TN_calls)
    FP <- sum(calls$rates$FP_calls)
    FN <- sum(calls$rates$FN_calls)
    
    bar_components <- rbind(bar_components,
                            data.frame(dataset = this_dataset, method = method, count = TP, type = "TP"))
    bar_components <- rbind(bar_components,
                            data.frame(dataset = this_dataset, method = method, count = FN, type = "FN"))
    bar_components <- rbind(bar_components,
                            data.frame(dataset = this_dataset, method = method, count = FP, type = "FP"))
    bar_components <- rbind(bar_components,
                            data.frame(dataset = this_dataset, method = method, count = TN, type = "TN"))
  }
  
  bar_components$dataset <- factor(bar_components$dataset, levels = datasets)
  levels(bar_components$dataset) <- paste0(levels(bar_components$dataset), " et al.")
  
  bar_components$type <- factor(bar_components$type, levels = c("TP", "TN", "FN", "FP"))
  # levels(bar_components$type) <- c("True positive", "True negative", "False negative", "False positive")
  
  p1 <- ggplot(bar_components, aes(x = method, y = count, fill = type, label = count)) +
    geom_bar(stat = "identity") +
    geom_text(size = 3, position = position_stack(vjust = 0.5), colour = "black", fontface = "bold") +
    theme_bw() +
    # scale_fill_brewer(palette = "GnBu") +
    scale_fill_manual(values = c("#ec561d", "#ffaf60", "#217db4", "#9bd4e4")) +
    labs(fill = "",
         x = "")
  
  # ----------------------------------------------------------------------------
  #   Plots 2-3: Parse predicted and observed outcomes
  # ----------------------------------------------------------------------------
  
  for(use_result_type in c("TPR", "FPR")) {
    legend <- NULL
    # Wrangle the poorly organized data
    plot_df <- results %>%
      filter(dataset == this_dataset) %>%
      filter(threshold == thresholds[[this_dataset]]) %>%
      filter(result_type == "predicted") %>%
      filter(score_type == use_result_type) %>%
      dplyr::select(!c(threshold, score_type))
    plot_df$true <- NA
    for(j in 1:nrow(plot_df)) {
      plot_df$true[j] <- results %>%
        filter(dataset == this_dataset) %>%
        filter(threshold == thresholds[[this_dataset]]) %>%
        filter(result_type == "true") %>%
        filter(score_type == use_result_type) %>%
        filter(DE_method == plot_df$DE_method[j]) %>%
        pull(point)
      
      p <- ggplot() +
        geom_boxplot(data = plot_df,
                     mapping = aes(x = DE_method,
                                   ymin = lower90,
                                   lower = lower50,
                                   middle = point,
                                   upper = upper50,
                                   ymax = upper90,
                                   fill = DE_method),
                     stat = "identity", color = "#666666", width = 0.5, alpha = 0.4) +
        geom_point(data = plot_df,
                   mapping = aes(x = DE_method, y = true, fill = DE_method),
                   size = 4, shape = 21, stroke = 1) +
        ylim(c(0,1)) +
        # scale_fill_manual(values = palette) +
        scale_fill_brewer(palette = "Set2") +
        theme_bw() +
        theme(legend.position = "none") +
        labs(x = "",
             y = ifelse(use_result_type == "TPR", "sensitivity", "specificity"))
      
      if(use_result_type == "TPR") {
        p2 <- p
      } else {
        p3 <- p
      }
    }
  }
  
  p2_padded <- plot_grid(p2, p3, ncol = 2, scale = 0.95)
  p <- plot_grid(p1, p2_padded, ncol = 2, rel_widths = c(1, 1.5))
  plotlist[[this_dataset]] <- p
}

common_scale <- 0.98

prow1 <- plot_grid(plotlist = plotlist[1:6],
                   ncol = 1,
                   labels = c("a", "b", "c", "d", "e", "f"),
                   scale = common_scale,
                   label_y = 1.04,
                   label_size = 18)

ggsave(file.path("output", "images", "results_summary_1.png"),
       plot = prow1,
       units = "in",
       height = 14,
       width = 9,
       bg = "white")

prow2 <- plot_grid(plotlist = plotlist[7:12],
                   ncol = 1,
                   labels = c("g", "h", "i", "j", "k", "m"),
                   scale = common_scale,
                   label_y = 1.04,
                   label_size = 18)

ggsave(file.path("output", "images", "results_summary_2.png"),
       plot = prow2,
       units = "in",
       height = 14,
       width = 9,
       bg = "white")
