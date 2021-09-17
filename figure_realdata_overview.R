library(tidyverse)
library(codaDE)
library(rulesoflife) # for generate_highcontrast_palette()

source("ggplot_fix.R")

palette_gen <- colorRampPalette(brewer.pal(9, "Set1")[1:8])
palette <- sample(palette_gen(100000))

dir.create("output", showWarnings = FALSE)
dir.create(file.path("output", "images"), showWarnings = FALSE)

use_sample_size <- 6

visualize_totals <- function(dataset, group_labels = NULL, min_relab = 0.01) {
  data <- do.call(paste0("parse_", dataset), list(absolute = TRUE))
  counts <- data$counts
  if(is.null(group_labels)) {
    groups <- factor(data$groups)
  } else {
    groups <- factor(data$groups, levels = names(group_labels))
    levels(groups) <- unname(unlist(group_labels))
  }
  
  g1_idx <- which(groups == unique(groups)[1])
  g2_idx <- which(groups == unique(groups)[2])
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

  counts <- counts[,idx_subset]
  groups <- groups[idx_subset]

  # Alternative version: select for plotting only the largest k features
  # k <- min(round(nrow(counts)*0.8), max_feature_n)  
  # feat_idx <- suppressMessages(data.frame(x = 1:nrow(counts), mu = rowMeans(counts)) %>%
  #                                arrange(desc(mu)) %>%
  #                                top_n(k) %>%
  #                                arrange(x) %>%
  #                                pull(x))
  # counts <- counts[feat_idx,]
  # count_sums <- colSums(counts[-feat_idx,])
  # counts <- rbind(counts, count_sums) # features x samples

  plots <- list()
  for(i in 1:2) {
    if(i == 1) {
      relab <- counts
      for(j in 1:ncol(counts)) {
        relab[,j] <- relab[,j] / sum(relab[,j])
      }
      mean_feature_relab <- rowMeans(relab)
      feature_palette <- sample(palette[1:nrow(relab)])
      filt <- mean_feature_relab < min_relab
      grays <- c(#"#808080",
                 #"#8a8a8a",
                 #"#949494",
                 # "#9d9d9d",
                 # "#a7a7a7",
                 # "#b1b1b1",
                 "#bbbbbb",
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
                              paste0(as.character(unique(groups)[1]), " (sample ", 1:length(g1_idx), ")"),
                              paste0(as.character(unique(groups)[2]), " (sample ", 1:length(g2_idx), ")"))
    
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
    plots[[i]] <- pl
  }
  pl <- plot_grid(plotlist = plots, ncol = 2) #, labels = c(alpha_label))
                     # #hjust = -1.35,
                     # vjust = -0.25, label_size = 24)
  return(pl)
}

plots <- list()
plots[[1]] <- visualize_totals(dataset = "Barlow",
                               group_labels = list(control = "Control diet", keto = "Ketogenic diet"))
plots[[2]] <- visualize_totals(dataset = "VieiraSilva",
                               group_labels = list(CD = "Crohn's disease", mHC = "Control"))
plots[[3]] <- visualize_totals(dataset = "Song",
                               group_labels = list(brain = "Brain", lung = "Lung"))
plots[[4]] <- visualize_totals(dataset = "Monaco",
                               group_labels = list(CD4_naive = "Naive CD4 cells", PBMC = "PBMC cells"))
plots[[5]] <- visualize_totals(dataset = "Hagai",
                               group_labels = list(pIC4 = "pIC4", unstimulated = "Unstimulated fibroblasts"))
plots[[6]] <- visualize_totals(dataset = "Owens",
                               group_labels = list(early_series = "Early time series", late_series = "Late time series"))
plots[[7]] <- visualize_totals(dataset = "Klein",
                               group_labels = list("LIF-2hr" = "LIF-treated", unstimulated = "Untreated"))
plots[[8]] <- visualize_totals(dataset = "Yu",
                               group_labels = list(Brn = "Brain", Lvr = "Liver"))

# Put these together
prow <- plot_grid(plotlist = plots, ncol = 1, labels = c("a", "b", "c", "d", "e", "f", "g", "h"))
ggsave(file.path("output", "images", "realdata_summary.png"),
       plot = prow,
       units = "in",
       height = 12,
       width = 8,
       bg = "white")
