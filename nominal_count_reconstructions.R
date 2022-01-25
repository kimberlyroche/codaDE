library(codaDE)
library(tidyverse)
library(cowplot)

# datasets <- c("VieiraSilva", "Barlow", "Song", "Monaco", "Hagai", "Owens", "Klein", "Yu")
# thresholds <- c(1, 1, 1, 1, 1, 1, 1, 1)

datasets <- c("Owens", "Klein", "Yu")
thresholds <- c(1, 1, 1)

plots <- list()
legends <- list()
for(i in 1:length(datasets)) {
  all_data <- readRDS(file.path("output", paste0("filtered_data_", datasets[i], "_threshold", thresholds[i], ".rds")))
  scaled_counts_A <- scaled_counts_ALDEx2(all_data$relative, pseudocount = 0.5)
  scaled_counts_D <- scaled_counts_DESeq2(all_data$relative, all_data$groups, pseudocount = 0.5)
  scaled_counts_S <- scaled_counts_scran(all_data$relative, all_data$groups, pseudocount = 0.5)
  
  plot_df <- data.frame(sample = 1:nrow(all_data$relative),
                        total = rowSums(all_data$relative),
                        group = all_data$groups,
                        type = "relative")
  plot_df <- rbind(plot_df,
                   data.frame(sample = 1:nrow(scaled_counts_A),
                              total = rowSums(scaled_counts_A),
                              group = all_data$groups,
                              type = "ALDEx2"))
  plot_df <- rbind(plot_df,
                   data.frame(sample = 1:nrow(scaled_counts_D),
                              total = rowSums(scaled_counts_D),
                              group = all_data$groups,
                              type = "DESeq2"))
  plot_df <- rbind(plot_df,
                   data.frame(sample = 1:nrow(scaled_counts_S),
                              total = rowSums(scaled_counts_S),
                              group = all_data$groups,
                              type = "scran"))
  plot_df <- rbind(plot_df,
                   data.frame(sample = 1:nrow(all_data$absolute),
                              total = rowSums(all_data$absolute),
                              group = all_data$groups,
                              type = "nominal"))
  
  plot_df$type <- factor(plot_df$type, levels = c("relative",
                                                  "nominal",
                                                  "ALDEx2",
                                                  "DESeq2",
                                                  "scran"))
  
  if(datasets[i] == "Owens") {
    plot_df$group <- factor(plot_df$group, levels = c("early_series", "late_series"))
    levels(plot_df$group) <- c("Early series", "Late series")
  } else if(datasets[i] == "Klein") {
    plot_df$group <- factor(plot_df$group, levels = c("unstimulated", "LIF-2hr"))
    levels(plot_df$group) <- c("Unstimulated", "LIF 2-hr treatment")
  } else if(datasets[i] == "Yu") {
    plot_df$group <- factor(plot_df$group, levels = c("Brn", "Lvr"))
    levels(plot_df$group) <- c("Brain", "Liver")
  }
  
  p <- ggplot(plot_df, aes(x = sample, y = total, fill = group)) +
    geom_bar(stat = "identity", alpha = 1) +
    facet_wrap(. ~ type, scales = "free_y", ncol = 5) +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    labs(fill = "Condition",
         x = "sample index",
         y = "total abundance")
  legend <- get_legend(p)
  p <- p +
    theme(legend.position = "none")
  
  plots[[i]] <- p
  legends[[i]] <- legend
}

l1 <- plot_grid(legends[[1]], NULL, ncol = 2, rel_widths = c(1, 0.4))
l2 <- plot_grid(legends[[2]], NULL, ncol = 2, rel_widths = c(1, 0.05))
l3 <- plot_grid(legends[[3]], NULL, ncol = 2, rel_widths = c(1, 0.8))
p1 <- plot_grid(plots[[1]], l1, ncol = 2, rel_widths = c(1, 0.17))
p2 <- plot_grid(plots[[2]], l2, ncol = 2, rel_widths = c(1, 0.17))
p3 <- plot_grid(plots[[3]], l3, ncol = 2, rel_widths = c(1, 0.17))
p <- plot_grid(p1, p2, p3, ncol = 1)

ggsave(file.path("output",
                 "images",
                 "nominal_rescaled.png"),
       p,
       dpi = 100,
       units = "in",
       height = 6,
       width = 12)
