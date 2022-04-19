source("path_fix.R")

library(tidyverse)
library(codaDE)
library(RSQLite)
library(gridExtra)
library(RColorBrewer)
library(cowplot)
library(ggridges)

palette_gen <- colorRampPalette(brewer.pal(9, "Set1")[1:8])
palette <- sample(palette_gen(100000))

use_sample_size <- 10

visualize_totals <- function(dataset, group_labels = NULL, min_relab = 0.01) {
  data <- readRDS(paste0("output/filtered_data_", dataset, "_threshold1.rds"))
  counts <- t(data$absolute)
  groups <- data$groups
  
  g1_idx <- which(groups == names(group_labels)[1])
  g2_idx <- which(groups == names(group_labels)[2])
  if(min(length(g1_idx), length(g2_idx)) < use_sample_size) {
    # Downsample to the smaller cohort size
    use_sample_size <- min(length(g1_idx), length(g2_idx))
  }
  if(length(g1_idx) > use_sample_size) {
    g1_idx <- sample(g1_idx, size = use_sample_size)
  }
  if(length(g2_idx) > use_sample_size) {
    g2_idx <- sample(g2_idx, size = use_sample_size)
  }
  
  idx_subset <- c(g1_idx, g2_idx)
  
  relab <- t(data$relative)
  for(j in 1:ncol(relab)) {
    relab[,j] <- relab[,j] / sum(relab[,j])
  }
  
  counts <- counts[,idx_subset]
  relab <- relab[,idx_subset]
  groups <- groups[idx_subset]
  
  plots <- list()
  for(i in 1:2) {
    if(i == 1) {
      mean_feature_relab <- rowMeans(relab)
      feature_palette <- sample(palette[1:nrow(relab)])
      filt <- mean_feature_relab < min_relab
      grays <- c("#bbbbbb",
                 "#c5c5c5",
                 "#cecece",
                 "#d8d8d8",
                 "#e2e2e2")
      feature_palette[filt] <- sample(grays, size = sum(filt), replace = TRUE)
    } else {
      # Relative abundances should already have been computed
      counts <- relab
    }
    
    n_samples <- ncol(counts)
    counts_idx <- cbind(1:nrow(counts), counts)
    colnames(counts_idx) <- c("feature",
                              paste0(as.character(names(group_labels)[1]), " (sample ", 1:length(g1_idx), ")"),
                              paste0(as.character(names(group_labels)[2]), " (sample ", 1:length(g2_idx), ")"))
    
    data_long <- pivot_longer(as.data.frame(counts_idx),
                              !feature,
                              names_to = "sample",
                              values_to = "abundance")
    data_long$feature <- factor(data_long$feature)
    data_long$sample <- data_long$sample
    data_long$group <- sapply(data_long$sample, function(x) {
      str_split(x, " \\(")[[1]][1]
    })
    data_long$sample <- unname(sapply(data_long$sample, function(x) {
      pieces <- str_match(x, "(.*?) \\((.*?)\\)")
      pieces[1,ncol(pieces)]
    }))
    data_long$sample <- factor(data_long$sample,
                               levels = paste0("sample ", 1:max(length(g1_idx), length(g2_idx))))
    
    tick_labels <- as.character(groups)
    names(tick_labels) <- colnames(counts_idx)[2:ncol(counts_idx)]
    
    data_long$group <- factor(data_long$group, levels = names(group_labels))
    levels(data_long$group) <- unname(unlist(group_labels))
    
    pl <- ggplot(data_long, aes(fill = feature, y = abundance, x = sample)) + 
      geom_bar(position = "stack", stat = "identity") +
      scale_fill_manual(values = feature_palette) +
      scale_x_discrete(labels = tick_labels) +
      facet_wrap(. ~ group) +
      theme_bw() +
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0))) +
      labs(y = ifelse(i == 1, "abundance", "relative abundance"))
    if(dataset == "Kimmerling" & i == 1) {
      # Force scientific notation
      pl <- pl +
        scale_y_continuous(breaks = c(0, 1000000, 2000000),
                           labels = c("0", "1e+06", "2e+06"))
    }
    plots[[i]] <- pl
  }
  pl <- plot_grid(plotlist = plots, ncol = 2) #, labels = c(alpha_label))
  return(pl)
}

# Note: these take at least 10 min. render
plots <- list()
plots[[1]] <- visualize_totals(dataset = "VieiraSilva",
                               group_labels = list(mHC = "Control", CD = "Crohn's disease"))
plots[[2]] <- visualize_totals(dataset = "Hagai",
                               group_labels = list(unstimulated = "Unstimulated fibroblasts", pIC4 = "pIC4"))
plots[[3]] <- visualize_totals(dataset = "Hashimshony",
                               group_labels = list("0" = "Quiescent", "1" = "Cycling"))
plots[[4]] <- visualize_totals(dataset = "Song",
                               group_labels = list(brain = "Brain metastasis", lung = "Lung primary tumor"))
plots[[5]] <- visualize_totals(dataset = "Monaco",
                               group_labels = list(CD4_naive = "Naive CD4 cells", PBMC = "PBMC cells"))
plots[[6]] <- visualize_totals(dataset = "Barlow",
                               group_labels = list(control = "Control diet", keto = "Ketogenic diet"))
plots[[7]] <- visualize_totals(dataset = "Gruen",
                               group_labels = list(A = "Two-inhibitor medium", B = "Serum"))
plots[[8]] <- visualize_totals(dataset = "Muraro",
                               group_labels = list(alpha = "Alpha", beta = "Beta"))
plots[[9]] <- visualize_totals(dataset = "Kimmerling",
                               group_labels = list(low_mass = "Low mass", high_mass = "High mass"))
plots[[10]] <- visualize_totals(dataset = "Yu",
                                group_labels = list(Brn = "Brain", Lvr = "Liver"))
plots[[11]] <- visualize_totals(dataset = "Owens",
                                group_labels = list(early_series = "Early time series", late_series = "Late time series"))
plots[[12]] <- visualize_totals(dataset = "Klein",
                                group_labels = list(unstimulated = "Untreated", "LIF-2hr" = "LIF-treated"))

# Put these together
# These are roughly ordered by percent of features differentially abundant
prow1 <- plot_grid(plotlist = plots[1:6], ncol = 1, labels = c("a", "b", "c", "d", "e", "f"))
prow2 <- plot_grid(plotlist = plots[7:12], ncol = 1, labels = c("g", "h", "i", "j", "k", "m"))

dev.off()
tiff(file.path("output", "images", "F4.tif"), units = "in", width = 8, height = 10, res = 300)
prow1
dev.off()

tiff(file.path("output", "images", "F5.tif"), units = "in", width = 8, height = 10, res = 300)
prow2
dev.off()
