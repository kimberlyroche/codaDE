source("path_fix.R")

library(tidyverse)
library(codaDE)
library(RSQLite)
library(gridExtra)
library(RColorBrewer)
library(cowplot)
library(rulesoflife) # for generate_highcontrast_palette()

source("ggplot_fix.R")

# Create directories manually

dir.create("output", showWarnings = FALSE)
dir.create(file.path("output", "images"), showWarnings = FALSE)

# ------------------------------------------------------------------------------
#
#   Supplemental figure 1: distributional characteristics of simulations
#
# ------------------------------------------------------------------------------

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

qq_res <- dbGetQuery(conn, paste0("SELECT FC_ABSOLUTE FROM datasets ",
                                  "WHERE FC_ABSOLUTE <= 10 AND ",
                                  "FC_ABSOLUTE >= 0.1"))$FC_ABSOLUTE
qq_res <- sapply(qq_res, function(x) {
  if(x < 1) {
    1 / x
  } else {
    x
  }
})
qq <- quantile(qq_res, probs = seq(from = 0, to = 1, length.out = 4))

ramped_palette <- colorRampPalette(c("#7479c4", "#e63030"))(3)

res <- dbGetQuery(conn, paste0("SELECT P, CORRP, PERCENT_DIFF_REALIZ, ",
                               "FC_ABSOLUTE, FC_PARTIAL FROM datasets ",
                               "WHERE FC_ABSOLUTE <= 10 AND FC_ABSOLUTE >= 0.1"))

dbDisconnect(conn)

temp <- res %>%
  filter(CORRP %in% c(0, 4))
temp$CORRP <- factor(temp$CORRP, levels = c("0", "4"))
levels(temp$CORRP) <- c("Independent features", "Strongly correlated features")

lr_margin <- 0.45

# ------------------------------------------------------------------------------
#  Percent differential features
# ------------------------------------------------------------------------------

p1 <- ggplot(temp, aes(x = PERCENT_DIFF_REALIZ, fill = factor(P))) +
  geom_density(alpha = 0.6) +
  scale_fill_brewer(palette = "RdYlGn") +
  # facet_wrap(. ~ CORRP) +
  xlim(c(0, 1)) +
  theme_bw() +
  labs(x = "percent differentially abundant features",
       fill = "Feature number") +
  theme(legend.position = "none",
        plot.margin = margin(t = 0.4, l = lr_margin, r = lr_margin, b = 0.5, "cm"))

# ------------------------------------------------------------------------------
#   Absolute fold change in totals
# ------------------------------------------------------------------------------

p2 <- ggplot(temp, aes(x = FC_ABSOLUTE, fill = factor(P))) +
  geom_density(alpha = 0.6) +
  scale_fill_brewer(palette = "RdYlGn") +
  # facet_wrap(. ~ CORRP) +
  xlim(c(0, 10)) +
  theme_bw() +
  labs(x = "fold change in abundance",
       fill = "Feature number") +
  theme(legend.position = "none",
        plot.margin = margin(t = 0.4, l = lr_margin, r = lr_margin, b = 0.5, "cm"))

# ------------------------------------------------------------------------------
#   Percent zeros
# ------------------------------------------------------------------------------

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
res <- dbGetQuery(conn, paste0("SELECT P, COMP_C_P0_A, COMP_C_P0_B ",
                               "FROM characteristics LEFT JOIN datasets ",
                               "ON characteristics.UUID=datasets.UUID ",
                               "WHERE PARTIAL = 0;"))
dbDisconnect(conn)

temp2 <- data.frame(sparsity = c(res$COMP_C_P0_A, res$COMP_C_P0_B),
                      P = c(res$P, res$P))
p3 <- ggplot(temp2, aes(x = sparsity, fill = factor(P))) +
  geom_density(alpha = 0.6) +
  scale_fill_brewer(palette = "RdYlGn") +
  xlim(c(0,1)) +
  theme_bw() +
  labs(x = "percent zeros",
       fill = "Feature number")
legend <- get_legend(p3 + theme(legend.position = "bottom"))
p3 <- p3  +
  theme(legend.position = "none",
        plot.margin = margin(t = 0.4, l = lr_margin, r = lr_margin, b = 0.5, "cm"))

prow <- plot_grid(p1,
                  p2,
                  p3,
                  align = 'vh',
                  labels = c("a", "b", "c"),
                  hjust = -1,
                  nrow = 1,
                  label_size = 18,
                  label_x = -0.03,
                  label_y = 1.03)
pl <- plot_grid(prow,
                legend,
                ncol = 1,
                rel_heights = c(1, 0.1))
show(pl)
ggsave(file.path("output",
                 "images",
                 paste0("simulation_characteristics.png")),
       plot = pl,
       dpi = 100,
       units = "in",
       height = 3,
       width = 9)

# ------------------------------------------------------------------------------
#
#   Supplemental figure 2: ALDEx2/DESeq2/scran sensitivity vs. NB GLM
#
# ------------------------------------------------------------------------------

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

res <- dbGetQuery(conn, paste0("SELECT ",
                               "datasets.UUID AS UUID, ",
                               "METHOD, ",
                               "PARTIAL_INFO, ",
                               "BASELINE_TYPE, ",
                               "datasets.BASELINE_CALLS AS ORACLE_BASELINE, ",
                               "CALLS, ",
                               "results.BASELINE_CALLS AS SELF_BASELINE, ",
                               "P, ",
                               "CORRP, ",
                               "LOG_MEAN, ",
                               "PERTURBATION, ",
                               "REP_NOISE, ",
                               "FC_ABSOLUTE, ",
                               "FC_RELATIVE, ",
                               "FC_PARTIAL, ",
                               "MED_ABS_TOTAL, ",
                               "MED_REL_TOTAL, ",
                               "PERCENT_DIFF_REALIZ, ",
                               "TPR, ",
                               "FPR ",
                               "FROM results LEFT JOIN datasets ON ",
                               "results.UUID=datasets.UUID ",
                               "WHERE PARTIAL_INFO=0 ",
                               "AND BASELINE_TYPE='self' ",
                               "AND FC_ABSOLUTE <= 10 ",
                               "AND FC_ABSOLUTE >= 0.1;"))

# Strip "result-less" entries
res <- res %>%
  filter(!is.na(TPR) & !is.na(FPR))

dbDisconnect(conn)

results <- NULL

for(method in c("ALDEx2", "DESeq2", "scran")) {
  cat(paste0("Making calls for method ", method, "\n"))
  # Filter to method of interest
  res_method <- res %>%
    filter(METHOD == method)
  # Pull oracle 'yes' total
  NB_calls <- unname(sapply(res_method$ORACLE_BASELINE,
                            function(x) {
                              sum(p.adjust(as.numeric(strsplit(x, ";")[[1]]), method = "BH") < 0.05)
                            }))
  # Pull method 'yes' total
  method_calls <- unname(sapply(res_method$SELF_BASELINE,
                                function(x) {
                                  sum(p.adjust(as.numeric(strsplit(x, ";")[[1]]), method = "BH") < 0.05)
                                }))
  results <- rbind(results,
                   data.frame(method = method,
                              oracle_total = NB_calls,
                              self_total = method_calls,
                              FPR = res_method$FPR))
}

# pl <- ggplot(results, aes(x = oracle_total, y = self_total, fill = FPR)) +
#   geom_point(size = 2, shape = 21) +
#   facet_wrap(. ~ method, ncol = 4) +
#   theme_bw() +
#   scale_fill_gradient2(low = "#bbbbbb", high = "red") +
#   labs(x = "number of differential features (NB GLM)",
#        y = "number of differential features (this method)") +
#   theme(axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
#         axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
#         plot.margin = margin(t = 0.25, l = 0.25, r = 0.25, b = 0.5, "cm"))
# show(pl)

p1 <- ggplot(results %>% filter(method == "ALDEx2"),
             aes(x = oracle_total, y = self_total, fill = FPR)) +
  geom_point(size = 2, shape = 21) +
  facet_wrap(. ~ method, ncol = 4) +
  theme_bw() +
  scale_fill_gradient2(low = "#bbbbbb", high = "red") +
  labs(x = "differential features (NB GLM)",
       y = "differential features (ALDEx2)") +
  theme(legend.position = "none") #+
  # theme(axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
  #       axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
  #       plot.margin = margin(t = 0.25, l = 0.25, r = 0.25, b = 0.5, "cm"))

p2 <- ggplot(results %>% filter(method == "DESeq2"),
             aes(x = oracle_total, y = self_total, fill = FPR)) +
  geom_point(size = 2, shape = 21) +
  facet_wrap(. ~ method, ncol = 4) +
  theme_bw() +
  scale_fill_gradient2(low = "#bbbbbb", high = "red") +
  labs(x = "differential features (NB GLM)",
       y = "differential features (DESeq2)") +
  theme(legend.position = "none") #+
# theme(axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
#       axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
#       plot.margin = margin(t = 0.25, l = 0.25, r = 0.25, b = 0.5, "cm"))

p3 <- ggplot(results %>% filter(method == "scran"),
             aes(x = oracle_total, y = self_total, fill = FPR)) +
  geom_point(size = 2, shape = 21) +
  facet_wrap(. ~ method, ncol = 4) +
  theme_bw() +
  scale_fill_gradient2(low = "#bbbbbb", high = "red") +
  labs(x = "differential features (NB GLM)",
       y = "differential features (scran)") #+
# theme(axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
#       axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
#       plot.margin = margin(t = 0.25, l = 0.25, r = 0.25, b = 0.5, "cm"))
legend <- get_legend(p3)
p3 <- p3 +
  theme(legend.position = "none") 

prow <- plot_grid(p1, p2, p3, legend, ncol = 4, rel_widths = c(1, 1, 1, 0.3))
prow

ggsave(file.path("output", "images", "methods_vs_NBGLM.png"),
       prow,
       dpi = 100,
       units = "in",
       height = 3,
       width = 10)

# Calculate R-squared values

for(this_method in c("ALDEx2", "DESeq2", "scran")) {
  cat(paste0("R^2 ", this_method, ": ",
             round(cor(results %>% filter(method == this_method) %>% pull(oracle_total),
                       results %>% filter(method == this_method) %>% pull(self_total))**2, 3), "\n"))
}

# ------------------------------------------------------------------------------
#
#   Supplemental figure 3: sensitivity x specificity plots for all methods vs.
#                          NB GLM calls on absolute abundance data
#
# ------------------------------------------------------------------------------

# See the file `figure_sim_results.R`

# ------------------------------------------------------------------------------
#
#   Supplemental figure 4: real data set visualization w/ stacked bar plots
#
# ------------------------------------------------------------------------------

palette_gen <- colorRampPalette(brewer.pal(9, "Set1")[1:8])
palette <- sample(palette_gen(100000))

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

