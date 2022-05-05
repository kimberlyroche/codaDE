source("path_fix.R")

library(codaDE)
library(ggplot2)
library(MASS)
library(cowplot)

plot_panel <- function(dataset_name, method, call, wrong_feature = NULL,
                       condition_labels = NULL) {
  # Parse calls
  # data <- readRDS(paste0("output/filtered_data_",dataset_name,"_threshold1.rds"))
  data <- wrangle_validation_data(dataset_name, threshold = 1, use_cpm = FALSE)
  calls_obj <- readRDS(paste0("output/real_data_calls/no_norm/calls_oracle_",method,"_", dataset_name, "_threshold1_noHKG.rds"))

  fgroups <- factor(data$groups)
  if(!is.null(condition_labels)) {
    levels(fgroups) <- condition_labels
  }
  palette <- c("#f25c4b", "#3199e8")
  names(palette) <- levels(fgroups)
  
  if(is.null(wrong_feature)) {
    if(call == "FN") {
      wrong_feature <- sample(which(calls_obj$rates$FN_calls), size = 1)
    } else {
      wrong_feature <- sample(which(calls_obj$rates$FP_calls), size = 1)
    }
  }
  A_idx <- sample(which(fgroups == levels(fgroups)[1]))
  B_idx <- sample(which(fgroups == levels(fgroups)[2]))
  
  # Plot nominal abundances
  plot_df <- data.frame(x = 1:length(A_idx),
                        y = data$ref_data[A_idx,wrong_feature],
                        condition = levels(fgroups)[1])
  plot_df <- rbind(plot_df,
                   data.frame(x = (length(A_idx)+1):(length(A_idx)+length(B_idx)),
                              y = data$ref_data[B_idx,wrong_feature],
                              condition = levels(fgroups)[2]))
  p1 <- ggplot(plot_df, aes(x = x, y = y, fill = condition)) +
    geom_point(size = 2, shape = 21) +
    scale_fill_manual(values = palette) +
    theme_bw() +
    labs(x = "sample index",
         y = "nominal abundance")
  
  # Plot relative abundances
  plot_df <- data.frame(x = 1:length(A_idx),
                        y = data$data[A_idx,wrong_feature],
                        condition = levels(fgroups)[1])
  plot_df <- rbind(plot_df,
                   data.frame(x = (length(A_idx)+1):(length(A_idx)+length(B_idx)),
                              y = data$data[B_idx,wrong_feature],
                              condition = levels(fgroups)[2]))
  p2 <- ggplot(plot_df, aes(x = x, y = y, fill = condition)) +
    geom_point(size = 2, shape = 21) +
    scale_fill_manual(values = palette) +
    theme_bw() +
    labs(x = "sample index",
         y = "relative abundance")
  
  # if(method == "DESeq2") {
  #   scaled_counts <- scaled_counts_DESeq2(data$relative, data$groups, 0.5)
  # } else if(method == "scran") {
  #   scaled_counts <- scaled_counts_scran(data$relative, data$groups, 0.5)
  # } else {
  #   scaled_counts <- scaled_counts_ALDEx2(data$relative, 0.5)
  # }
  
  # Plot method-rescaled abundances
  # plot_df <- data.frame(x = 1:length(A_idx),
  #                       y = scaled_counts[A_idx,wrong_feature],
  #                       condition = levels(fgroups)[1])
  # plot_df <- rbind(plot_df,
  #                  data.frame(x = (length(A_idx)+1):(length(A_idx)+length(B_idx)),
  #                             y = scaled_counts[B_idx,wrong_feature],
  #                             condition = levels(fgroups)[2]))
  # p3 <- ggplot(plot_df, aes(x = x, y = y, fill = condition)) +
  #   geom_point(size = 2, shape = 21) +
  #   scale_fill_manual(values = palette) +
  #   theme_bw() +
  #   labs(x = "sample index",
  #        y = paste0(method, "-scaled abundance"))
  
  legend <- get_legend(p1 + theme(legend.position = "bottom"))
  
  prow <- plot_grid(p1 + theme(legend.position = "none"),
                    p2 + theme(legend.position = "none"),
                    ncol = 2,
                    scale = 0.92)
  
  p <- plot_grid(prow,
                 legend,
                 ncol = 1,
                 rel_heights = c(1, 0.05),
                 # labels = c("a", "b", "c"),
                 # label_size = 18,
                 scale = 0.92)
  p
}

panel_Kimmerling <- plot_panel(dataset_name = "Kimmerling", method = "DESeq2", call = "FN", wrong_feature = 3300,
                               condition_labels = c("low mass cells", "high mass cells")) # gene Ung

panel_Klein <- plot_panel(dataset_name = "Klein", method = "scran", call = "FP", wrong_feature = 1076) # gene H2-T3

p <- plot_grid(panel_Kimmerling, panel_Klein, ncol = 1,
               labels = c("a", "b"),
               label_size = 18,
               label_y = 1.03)

dev.off()
tiff(file.path("output", "images", "sample_calls.tif"), units = "in", width = 6, height = 6, res = 300)
p
dev.off()

# ggsave(file.path("output", "images", "spurious_calls.svg"),
#        p,
#        dpi = 100,
#        units = "in",
#        height = 6,
#        width = 9)

# Re-run the test to make sure this is "not differential" by the NB GLM
# call_DA_NB(data = data$absolute[,wrong_feature,drop=FALSE], fgroups)
