source("path_fix.R")

library(tidyverse)
library(codaDE)
library(cowplot)
library(RColorBrewer)

# Create directories manually
dir.create("output", showWarnings = FALSE)
dir.create(file.path("output", "images"), showWarnings = FALSE)

# "Canonical" ordering
datasets <- c("VieiraSilva", "Hagai", "Hashimshony", "Song", "Monaco", "Barlow",
              "Gruen", "Muraro", "Kimmerling", "Yu", "Owens", "Klein")
thresholds <- c(1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1)

names(thresholds) <- datasets
model_dir <- "oracle"
# submodel_dir <- "regression_cpm"
submodel_dir <- "regression_plain"

palette <- list(ALDEx2 = "#46A06B",
                DESeq2 = "#FF5733",
                edgeR_TMM = "#7E54DE",
                'edgeR (TMM)' = "#7E54DE",
                scran = "#E3C012")

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

# results$DE_method <- factor(results$DE_method)
# levels(results$DE_method)[3] <- "edgeR (TMM)"

# R^2 for prediction
temp <- results %>% filter(score_type == "TPR")
cat(paste0("R^2 sensitivity: ", round(cor(temp %>% filter(result_type == "true") %>% pull(point),
                                          temp %>% filter(result_type == "predicted") %>% pull(point))**2, 3), "\n"))
temp <- results %>% filter(score_type == "FPR")
cat(paste0("R^2 specificity: ", round(cor(temp %>% filter(result_type == "true") %>% pull(point),
                                          temp %>% filter(result_type == "predicted") %>% pull(point))**2, 3), "\n"))

# Statistics related to upper bound on sensitivity
temp <- results %>% filter(score_type == "TPR")
temp_true <- temp %>%
  filter(result_type == "true") %>%
  select(dataset, DE_method, point)
temp_pred <- temp %>%
  filter(result_type == "predicted") %>%
  select(dataset, DE_method, lower90)
temp2 <- temp_pred %>%
  left_join(temp_true, by = c("dataset", "DE_method")) %>%
  mutate(bounded = (point > lower90))
sum(temp2$bounded)/nrow(temp2)

# Statistics related to lower bound on specificity
temp <- results %>% filter(score_type == "FPR")
temp_true <- temp %>%
  filter(result_type == "true") %>%
  select(dataset, DE_method, point)
temp_pred <- temp %>%
  filter(result_type == "predicted") %>%
  select(dataset, DE_method, lower90)
temp2 <- temp_pred %>%
  left_join(temp_true, by = c("dataset", "DE_method")) %>%
  mutate(bounded = (point > lower90))
sum(temp2$bounded)/nrow(temp2)

plotlist <- list()
for(i in 1:length(datasets)) {
  this_dataset <- datasets[i]
  
  # ----------------------------------------------------------------------------
  #   Plot 1: Parse absolute counts of true positives, false positives, etc.
  # ----------------------------------------------------------------------------
  
  bar_components <- NULL
  for(method in c("ALDEx2", "DESeq2", "edgeR_TMM", "scran")) {
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
  
  bar_components$method <- factor(bar_components$method)
  levels(bar_components$method)[3] <- "edgeR (TMM)"
  
  p1 <- ggplot(bar_components, aes(x = method, y = count, fill = type, label = count)) +
    geom_bar(stat = "identity") +
    geom_text(size = 3, position = position_stack(vjust = 0.5), colour = "black", fontface = "bold") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
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
    }
      
    plot_df$DE_method <- factor(plot_df$DE_method)
    levels(plot_df$DE_method)[3] <- "edgeR (TMM)"
    
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
                 size = 3, shape = 21, stroke = 1) +
      ylim(c(0,1)) +
      scale_fill_manual(values = palette) +
      # scale_fill_brewer(palette = "Set2") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      labs(x = "",
           y = ifelse(use_result_type == "TPR", "sensitivity", "specificity"))
    
    if(use_result_type == "TPR") {
      p2 <- p
    } else {
      p3 <- p
    }
  }
  
  p2_padded <- plot_grid(p2, p3, ncol = 2, scale = 0.95)
  p <- plot_grid(p1, p2_padded, ncol = 2, rel_widths = c(1, 1.5))
  plotlist[[this_dataset]] <- p
}

common_scale <- 0.98

neg_pad <- -0.12

prow1 <- plot_grid(plotlist[[1]], NULL,
                   plotlist[[2]], NULL,
                   plotlist[[3]], NULL,
                   plotlist[[4]], NULL,
                   plotlist[[5]], NULL,
                   plotlist[[6]],
                   ncol = 1,
                   rel_heights = c(1, neg_pad, 1, neg_pad, 1, neg_pad,
                                   1, neg_pad, 1, neg_pad, 1),
                   labels = c("a", "", "b", "", "c", "", "d", "", "e", "", "f"),
                   scale = common_scale,
                   label_y = 1.04,
                   label_size = 18)

tiff(file.path("output", "images", "F6.tif"), units = "in", width = 9, height = 14, res = 300)
prow1
dev.off()

prow2 <- plot_grid(plotlist[[7]], NULL,
                   plotlist[[8]], NULL,
                   plotlist[[9]], NULL,
                   plotlist[[10]], NULL,
                   plotlist[[11]], NULL,
                   plotlist[[12]],
                   ncol = 1,
                   rel_heights = c(1, neg_pad, 1, neg_pad, 1, neg_pad,
                                   1, neg_pad, 1, neg_pad, 1),
                   labels = c("g", "", "h", "", "i", "", "j", "", "k", "", "l"),
                   scale = common_scale,
                   label_y = 1.04,
                   label_size = 18)


tiff(file.path("output", "images", "F7.tif"), units = "in", width = 9, height = 14, res = 300)
prow2
dev.off()

