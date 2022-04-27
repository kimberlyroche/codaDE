source("path_fix.R")

library(tidyverse)
library(codaDE)
library(RSQLite)
library(cowplot)
library(ggExtra)
library(RColorBrewer)
library(randomForest)
library(ggridges)

output_dir <- file.path("output", "images")

# This is only in use in some commented out code below
string_to_calls <- function(call_string) {
  pvals <- as.numeric(strsplit(call_string, ";")[[1]])
  pvals <- p.adjust(pvals, method = "BH")
  pvals < 0.05
}

pull_data <- function(qq, use_baseline = "self", use_cpm = TRUE) {
  conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
  res <- dbGetQuery(conn, paste0("SELECT ",
                                 "datasets.UUID AS UUID, ",
                                 "METHOD, ",
                                 "PARTIAL_INFO, ",
                                 "BASELINE_TYPE, ",
                                 "datasets.BASELINE_CALLS AS ORACLE_BASELINE, ",
                                 "datasets.BASELINE_BETAS AS ORACLE_BETAS, ",
                                 "CALLS, ",
                                 "results.BASELINE_CALLS AS SELF_BASELINE, ",
                                 "BETAS, ",
                                 "results.BASELINE_BETAS AS SELF_BETAS, ",
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
                                 "FPR, ",
                                 "FDR, ",
                                 "BETA ",
                                 "FROM results LEFT JOIN datasets ON ",
                                 "results.UUID=datasets.UUID ",
                                 "WHERE PARTIAL_INFO=0 ",
                                 "AND BASELINE_TYPE='", use_baseline, "' ",
                                 ifelse(use_cpm, "AND OBSERVED_TYPE='cpm' ", "AND OBSERVED_TYPE != 'cpm' "),
                                 "AND FC_ABSOLUTE <= 10 ",
                                 "AND FC_ABSOLUTE >= 0.1;"))
  dbDisconnect(conn)
  
  # Strip "result-less" entries. These are simulations with negligible amounts
  # of differential abundance - the simulated effect sizes were just too small!
  res <- res %>%
    filter(!is.na(TPR) & !is.na(FPR))
  
  # Filter out MAST entries. The added zero-inflation component makes MAST really
  # different than the other models and exceptionally sensitive to data processing
  # (esp filtering) choices. Remove for now.
  res <- res %>%
    filter(METHOD != "MAST")
  
  # Generate an absolute fold change low/med/high factor labeling
  # This reports the scale of increases in fold abundance either from A -> B or
  # from B -> A
  res$FC_plot <- sapply(res$FC_ABSOLUTE, function(x) {
    if(x < 1) {
      1 / x
    } else {
      x
    }
  })
  res$FC_plot <- cut(res$FC_plot, breaks = qq)
  levels(res$FC_plot) <- c(paste0("low (< ", qq[2], ")"),
                           "moderate",
                           paste0("high (> ", qq[3], ")"),
                           paste0("very high (> ", qq[4], ")"))
  
  # Generate an fold change direction label
  res$FC_dir <- sapply(res$FC_ABSOLUTE, function(x) {
    if(x < 1) {
      -1
    } else {
      1
    }
  })
  res$FC_dir <- factor(res$FC_dir)
  res
}

plot_ROC_flag <- function(data, method, p, logical_vec) {
  data$flag <- FALSE
  data$flag[logical_vec] <- TRUE
  pl <- ggplot() +
    geom_point(data = data %>% filter(!flag & METHOD == method & P == p),
               mapping = aes_string(x = "FPR", y = "TPR"),
               # size = 2,
               # shape = 21,
               color = "#dddddd") +
    geom_point(data = data %>% filter(flag & METHOD == method & P == p),
               mapping = aes_string(x = "FPR", y = "TPR"),
               size = 1.5,
               shape = 21,
               fill = "#1ab079") +
    xlim(c(0,1)) +
    ylim(c(0,1)) +
    labs(x = "FPR",
         y = "TPR") +
    theme_bw() +
    facet_wrap(. ~ METHOD, ncol = 4)
  return(pl)
}

plot_ROC_fc <- function(data, method, p = NULL) {
  pl <- ggplot()
  if(is.null(p)) {
    pl <- pl +
      geom_point(data %>% filter(METHOD == method),
                 mapping = aes(x = specificity, y = sensitivity, fill = FC_plot),
                 size = 1.5,
                 shape = 21)
  } else {
    pl <- pl +
      geom_point(data %>% filter(METHOD == method & P == p),
                 mapping = aes(x = specificity, y = sensitivity, fill = FC_plot),
                 size = 1.5,
                 shape = 21)
  }
  pl <- pl +
    xlim(c(-0.025,1.025)) +
    ylim(c(-0.025,1.025)) +
    labs(x = paste0(method, " specificity (1-FPR)"),
         y = paste0(method, " sensitivity"),
         fill = "Fold change") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 9, margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(size = 9, margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0))) +
    scale_fill_brewer(palette = "RdYlBu")
  legend <- get_legend(pl + theme(legend.position = "bottom"))
  pl <- pl +
    theme(legend.position = "none")
  pl <- ggMarginal(pl, margins = "both", type = "histogram", bins = 30,
                   fill = "#888888", color = "#333333", size = 10)
  return(list(pl = pl, legend = legend))
}

plot_ROC_percentDE <- function(data, method, p = NULL) {
  pl <- ggplot()
  if(is.null(p)) {
    pl <- pl +
      geom_point(data %>% filter(METHOD == method),
                 mapping = aes(x = specificity, y = sensitivity, fill = PERCENT_DIFF_REALIZ),
                 size = 1.5,
                 shape = 21)
  } else {
    pl <- pl +
      geom_point(data %>% filter(METHOD == method & P == p),
                 mapping = aes(x = specificity, y = sensitivity, fill = PERCENT_DIFF_REALIZ),
                 size = 1.5,
                 shape = 21)
  }
  pl <- pl +
    xlim(c(-0.025,1.025)) +
    ylim(c(-0.025,1.025)) +
    labs(x = paste0(method, " specificity (1-FPR)"),
         y = paste0(method, " sensitivity"),
         fill = "Proportion differentially abundant features     ") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 9, margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(size = 9, margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
          legend.position = "bottom",
          legend.title = element_text(vjust = 0.75)) +
    scale_fill_distiller(palette = "RdYlBu",
                         limits = c(0,1),
                         breaks = c(0, 0.5, 1),
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))
  legend <- get_legend(pl)
  pl <- pl +
    theme(legend.position = "none")
  pl <- ggMarginal(pl, margins = "both", type = "histogram", bins = 30,
                   fill = "#888888", color = "#333333", size = 10)
  return(list(pl = pl, legend = legend))
}

# Method-labeling palette
palette <- c("#46A06B", "#FF5733", "#EF82BB", "#9B59B6", "#E3C012", "#C70039", "#467BBC")
names(palette) <- c("ALDEx2", "ANCOM-BC", "DESeq2", "edgeR (TMM)", "scran", "DESeq2 (CG)", "edgeR")

# Visualize this palette
# ggplot(data.frame(x = 1:7, y = 1:7, method = names(palette)), aes(x = x, y = y, fill = method)) +
#   geom_point(size = 5, shape = 21) +
#   scale_fill_manual(values = palette) +
#   theme_bw()

# ------------------------------------------------------------------------------
#   Initialize output directory and generate a color palette for levels of fold
#   change
# ------------------------------------------------------------------------------

# Create directories manually if they don't already exist
dir.create("output", showWarnings = FALSE)
dir.create(file.path("output", "images"), showWarnings = FALSE)

method_list <- c("ALDEx2", "DESeq2", "scran", "edgeR", "edgeR_TMM", "ANCOMBC", "DESeq2_control")

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

qq <- c(0, 1.5, 2.5, 5, Inf)

# ------------------------------------------------------------------------------
#   Sensitivity x specificity plot labeled by fold change
#
#   Sensitivity and specificity of calls made by all methods on absolute vs.
#   relative abundances
#
#   For some reason, this code doesn't work as a function (I guess something 
#   needs to be global?) so it's duplicated below.
# ------------------------------------------------------------------------------

# data <- pull_data(qq, use_baseline = "oracle", use_cpm = TRUE)
data <- pull_data(qq, use_baseline = "oracle", use_cpm = FALSE)

data$sensitivity <- data$TPR
data$specificity <- 1 - data$FPR

data$METHOD <- factor(data$METHOD)
levels(data$METHOD)[2] <- "ANCOM-BC"
levels(data$METHOD)[4] <- "DESeq2 (CG)"
levels(data$METHOD)[6] <- "edgeR (TMM)"

# ------------------------------------------------------------------------------
#   Calculate some overall statistics
# ------------------------------------------------------------------------------

if(FALSE) {
  # k
  threshold <- 0.5
  
  # Exceeding k% FPR by SETTING
  data %>%
    select(METHOD, P, PERTURBATION, FPR) %>%
    mutate(setting = ifelse(P <= 1000 & PERTURBATION >= 1.2,
                            "1: Microbial",
                            ifelse(P >= 1000 & PERTURBATION >= 0.8 & PERTURBATION <= 1.2,
                                   "2: Cell transcriptomic",
                                   ifelse(P == 5000 & PERTURBATION <= 0.8,
                                          "3: Bulk transcriptomic",
                                          NA)))) %>%
    mutate(high_FPR = ifelse(FPR > threshold, 1, 0)) %>%
    mutate(low_FPR = ifelse(FPR <= threshold, 1, 0)) %>%
    group_by(METHOD, setting) %>%
    summarize(median_specificity = 1 - median(FPR),
              n_high_FPR = sum(high_FPR),
              n_low_FPR = sum(low_FPR)) %>%
    mutate(prop = round(n_high_FPR / (n_high_FPR + n_low_FPR), 2)) %>%
    arrange(setting, METHOD) %>%
    filter(METHOD != "edgeR" & !is.na(setting))
  
  # Exceeding k% FPR by FEATURE NUMBER
  data %>%
    select(METHOD, P, FPR) %>%
    mutate(high_FPR = ifelse(FPR > threshold, 1, 0)) %>%
    mutate(low_FPR = ifelse(FPR <= threshold, 1, 0)) %>%
    group_by(METHOD, P) %>%
    summarize(median_specificity = 1 - median(FPR),
              n_high_FPR = sum(high_FPR),
              n_low_FPR = sum(low_FPR)) %>%
    mutate(prop = round(n_high_FPR / (n_high_FPR + n_low_FPR), 2)) %>%
    arrange(P, METHOD) %>%
    filter(METHOD != "edgeR")
  
  # # 80% of data sets with modest change had FPR of 10% or less
  # quantile(data %>%
  #            filter(PERCENT_DIFF_REALIZ < 1/3) %>%
  #            filter(FC_plot == "low (< 1.5)") %>%
  #            pull(FPR), probs = c(0.8))
  # 
  # # Median specificity per method
  # data %>%
  #   group_by(METHOD) %>%
  #   summarize("Median sensitivity" = median(TPR),
  #             "Median specificity" = 1 - median(FPR))
  # 
  # # Exceeding k% FPR overall
  # threshold <- 0.50
  # data %>%
  #   select(METHOD, FPR) %>%
  #   mutate(high_FPR = ifelse(FPR > threshold, 1, 0)) %>%
  #   mutate(low_FPR = ifelse(FPR <= threshold, 1, 0)) %>%
  #   group_by(METHOD) %>%
  #   summarize(n_high_FPR = sum(high_FPR),
  #             n_low_FPR = sum(low_FPR)) %>%
  #   mutate(prop = round(n_high_FPR / (n_high_FPR + n_low_FPR), 2))
}

# ------------------------------------------------------------------------------
#   ROC plot: edgeR without normalization is terrible
# ------------------------------------------------------------------------------

p <- ggplot(data %>% filter(METHOD == "DESeq2"),
            aes(x = FPR, y = TPR, fill = factor(P))) +
  geom_point(size = 2, shape = 21) +
  theme_bw() +
  labs(x = "false positive rate",
       y = "true positive rate",
       fill = "feature number") +
  scale_fill_brewer(palette = "Reds") +
  theme(legend.position = "bottom")
pm <- ggMarginal(p, margins = "both", type = "histogram", bins = 50,
                 fill = "#888888", color = "#333333", size = 20)

ggsave(file.path(output_dir, "roc_edgeR-only.svg"),
       pm,
       dpi = 100,
       units = "in",
       height = 6,
       width = 6)

# ------------------------------------------------------------------------------
#   ROC plot: all normalization-based methods together
# ------------------------------------------------------------------------------

p <- ggplot(data %>% filter(METHOD %in% method_list),
            aes(x = FPR, y = TPR, fill = factor(P))) +
  geom_point(size = 2, shape = 21) +
  theme_bw() +
  labs(x = "false positive rate",
       y = "true positive rate",
       fill = "feature number") +
  scale_fill_brewer(palette = "Reds") +
  theme(legend.position = "bottom")
pm <- ggMarginal(p, margins = "both", type = "histogram", bins = 50,
                 fill = "#888888", color = "#333333", size = 20)

ggsave(file.path(output_dir, "roc_edgeR-absent.svg"),
       pm,
       dpi = 100,
       units = "in",
       height = 6,
       width = 6)

# ------------------------------------------------------------------------------
#   Ridges plot: all normalization-based methods x FPR
# ------------------------------------------------------------------------------

# Main text version
use_fdr <- 0.01

# Supplemental, higher-stringency version
# use_fdr <- 0.001

method_list <- c("ALDEx2", "ANCOM-BC", "DESeq2", "edgeR (TMM)", "scran")
temp <- data %>%
  filter(METHOD %in% method_list) %>%
  filter(FDR == use_fdr)
temp$METHOD <- factor(temp$METHOD, levels = method_list[5:1])

p <- ggplot(temp, aes(x = FPR, y = METHOD, fill = METHOD)) +
  geom_density_ridges(stat = "binline", bins = 40, scale = 1) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  scale_y_discrete(expand = expansion(add = c(0.2, 1.1))) +
  coord_cartesian(clip = "off") +
  guides(fill = "none") +
  labs(x = "false positive rate") +
  scale_fill_manual(values = palette)

ggsave(file.path(output_dir, "ridges_fpr_methods.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 4,
       width = 4.5)

# ------------------------------------------------------------------------------
#   Ridges plot: FPR x P
# ------------------------------------------------------------------------------

p <- ggplot(temp, aes(x = FPR, y = factor(P))) +
  geom_density_ridges(stat = "binline", bins = 40, scale = 1) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  scale_y_discrete(expand = expansion(add = c(0.2, 1.1))) +
  coord_cartesian(clip = "off") +
  guides(fill = "none") +
  labs(x = "false positive rate")

ggsave(file.path(output_dir, "ridges_fpr_P.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 4,
       width = 4)

# ------------------------------------------------------------------------------
#   Sensitivity x specificity plot labeled by percent of features differentially
#   abundant
# ------------------------------------------------------------------------------

by_pde <- TRUE

p1_100 <- NULL; p1_1000 <- NULL; p1_5000 <- NULL
p2_100 <- NULL; p2_1000 <- NULL; p2_5000 <- NULL
p3_100 <- NULL; p3_1000 <- NULL; p3_5000 <- NULL
p4_100 <- NULL; p4_1000 <- NULL; p4_5000 <- NULL
p5_100 <- NULL; p5_1000 <- NULL; p5_5000 <- NULL
legend <- NULL
for(i in 1:length(method_list)) {
  method <- method_list[i]
  for(j in c(100, 1000, 5000)) {
    if(by_pde) {
      plot_pieces <- plot_ROC_percentDE(temp, method, j)
    } else {
      plot_pieces <- plot_ROC_fc(temp, method, j) # fold change version
    }
    pl <- plot_pieces$pl
    if(is.null(legend)) {
      legend <- plot_pieces$legend
    }
    if(i == 1 & j == 100) p1_100 <<- pl
    if(i == 1 & j == 1000) p1_1000 <<- pl
    if(i == 1 & j == 5000) p1_5000 <<- pl
    if(i == 2 & j == 100) p2_100 <<- pl
    if(i == 2 & j == 1000) p2_1000 <<- pl
    if(i == 2 & j == 5000) p2_5000 <<- pl
    if(i == 3 & j == 100) p3_100 <<- pl
    if(i == 3 & j == 1000) p3_1000 <<- pl
    if(i == 3 & j == 5000) p3_5000 <<- pl
    if(i == 4 & j == 100) p4_100 <<- pl
    if(i == 4 & j == 1000) p4_1000 <<- pl
    if(i == 4 & j == 5000) p4_5000 <<- pl
    if(i == 5 & j == 100) p5_100 <<- pl
    if(i == 5 & j == 1000) p5_1000 <<- pl
    if(i == 5 & j == 5000) p5_5000 <<- pl
  }
}
prow1 <- plot_grid(p1_100, p1_1000, p1_5000, ncol = 3)
prow2 <- plot_grid(p2_100, p2_1000, p2_5000, ncol = 3)
prow3 <- plot_grid(p3_100, p3_1000, p3_5000, ncol = 3)
prow4 <- plot_grid(p4_100, p4_1000, p4_5000, ncol = 3)
prow5 <- plot_grid(p5_100, p5_1000, p5_5000, ncol = 3)

pgrid <- plot_grid(prow1, prow2, prow3, prow4, prow5,
                   nrow = 5, labels = c("a", "b", "c", "d", "e"),
                   label_size = 17, label_y = 1.02)
pl <- plot_grid(pgrid, legend, ncol = 1, rel_heights = c(1, 0.08))

fn <- ""
if(by_pde & use_fdr == 0.05) {
  fn <- "S3.svg"
} else if(!by_pde & use_fdr == 0.05) {
  fn <- "S4.svg"
} else if(by_pde) {
  fn <- "S3_alt-fdr.svg"
} else {
  fn <- "S4_alt-fdr.svg"
}
ggsave(file.path("output", "images", fn),
       plot = pl,
       units = "in",
       height = 12,
       width = 7.75,
       bg = "white")

# ------------------------------------------------------------------------------
#   Same plot by "settings"
# ------------------------------------------------------------------------------

by_pde <- TRUE

p1_1 <- NULL; p1_2 <- NULL; p1_3 <- NULL
p2_1 <- NULL; p2_2 <- NULL; p2_3 <- NULL
p3_1 <- NULL; p3_2 <- NULL; p3_3 <- NULL
p4_1 <- NULL; p4_2 <- NULL; p4_3 <- NULL
p5_1 <- NULL; p5_2 <- NULL; p5_3 <- NULL
legend <- NULL
for(i in 1:length(method_list)) {
  method <- method_list[i]
  setting_label <- "n/a"
  for(j in 1:3) { # iterate settings
    subdata <- temp
    if(j == 1) {
      subdata <- subdata %>%
        filter(P <= 1000 & PERTURBATION >= 1.2)
      setting_label <- "Microbial"
    } else if(j == 2) {
      subdata <- subdata %>%
        filter(P >= 1000 & PERTURBATION >= 0.8 & PERTURBATION <= 1.2)
      setting_label <- "Transcriptomic"
    } else if(j == 3) {
      subdata <- subdata %>%
        filter(P == 5000 & PERTURBATION <= 0.8)
      setting_label <- "Bulk transcriptomic"
    }
    prop_high_fpr <- subdata %>%
      filter(METHOD == method) %>%
      mutate(high_FPR = ifelse(FPR < 0.95, "high", "low")) %>%
      group_by(high_FPR) %>%
      tally() %>%
      mutate(prop = n/sum(n)) %>%
      filter(high_FPR == "high") %>%
      pull(prop)
    if(by_pde) {
      plot_pieces <- plot_ROC_percentDE(subdata, method)
    } else {
      plot_pieces <- plot_ROC_fc(subdata, method)
    }
    pl <- plot_pieces$pl
    if(is.null(legend)) {
      legend <- plot_pieces$legend
    }
    if(i == 1 & j == 1) p1_1 <<- pl
    if(i == 1 & j == 2) p1_2 <<- pl
    if(i == 1 & j == 3) p1_3 <<- pl
    if(i == 2 & j == 1) p2_1 <<- pl
    if(i == 2 & j == 2) p2_2 <<- pl
    if(i == 2 & j == 3) p2_3 <<- pl
    if(i == 3 & j == 1) p3_1 <<- pl
    if(i == 3 & j == 2) p3_2 <<- pl
    if(i == 3 & j == 3) p3_3 <<- pl
    if(i == 4 & j == 1) p4_1 <<- pl
    if(i == 4 & j == 2) p4_2 <<- pl
    if(i == 4 & j == 3) p4_3 <<- pl
    if(i == 5 & j == 1) p5_1 <<- pl
    if(i == 5 & j == 2) p5_2 <<- pl
    if(i == 5 & j == 3) p5_3 <<- pl
  }
}

prow1 <- plot_grid(p1_1, p1_2, p1_3, ncol = 3)
prow2 <- plot_grid(p2_1, p2_2, p2_3, ncol = 3)
prow3 <- plot_grid(p3_1, p3_2, p3_3, ncol = 3)
prow4 <- plot_grid(p4_1, p4_2, p4_3, ncol = 3)
prow5 <- plot_grid(p5_1, p5_2, p5_3, ncol = 3)
pgrid <- plot_grid(prow1, prow2, prow3, prow4, prow5,
                   nrow = 5, labels = c("a", "b", "c", "d", "e"),
                   label_size = 17, label_y = 1.02)
pl <- plot_grid(pgrid, legend, ncol = 1, rel_heights = c(1, 0.08))
ggsave(file.path("output", "images", ifelse(by_pde, "F1_alt.svg", "F2.svg")),
       plot = pl,
       units = "in",
       height = 12,
       width = 7.75,
       bg = "white")

# ------------------------------------------------------------------------------
#   DESeq2 improvement with control_genes (simulated data)
# ------------------------------------------------------------------------------

p1_1 <- NULL; p1_2 <- NULL; p1_3 <- NULL
p2_1 <- NULL; p2_2 <- NULL; p2_3 <- NULL
legend <- NULL
setting_label <- "n/a"
for(i in 1:2) {
  this_method <- c("DESeq2", "DESeq2 (CG)")[i]
  for(j in 1:3) { # iterate settings
    subdata <- data %>%
      filter(METHOD == this_method) %>%
      filter(FDR == 0.05)
    if(j == 1) {
      subdata <- subdata %>%
        filter(P <= 1000 & PERTURBATION >= 1.2)
      setting_label <- "Microbial"
    } else if(j == 2) {
      subdata <- subdata %>%
        filter(P >= 1000 & PERTURBATION >= 0.8 & PERTURBATION <= 1.2)
      setting_label <- "Transcriptomic"
    } else if(j == 3) {
      subdata <- subdata %>%
        filter(P == 5000 & PERTURBATION <= 0.8)
      setting_label <- "Bulk transcriptomic"
    }
    if(by_pde) {
      plot_pieces <- plot_ROC_percentDE(subdata, this_method)
    } else {
      plot_pieces <- plot_ROC_fc(subdata, this_method)
    }
    pl <- plot_pieces$pl
    if(is.null(legend)) {
      legend <- plot_pieces$legend
    }
    if(i == 1 & j == 1) p1_1 <<- pl
    if(i == 1 & j == 2) p1_2 <<- pl
    if(i == 1 & j == 3) p1_3 <<- pl
    if(i == 2 & j == 1) p2_1 <<- pl
    if(i == 2 & j == 2) p2_2 <<- pl
    if(i == 2 & j == 3) p2_3 <<- pl
  }
}

prow1 <- plot_grid(p1_1, p1_2, p1_3, ncol = 3)
prow2 <- plot_grid(p2_1, p2_2, p2_3, ncol = 3)

ptop <- plot_grid(prow1, prow2,
                  nrow = 2, labels = c("a", "b"),
                  label_size = 18, label_y = 1.02)

nocg_res <- data %>%
  filter(METHOD == "DESeq2" & FDR == 0.05) %>%
  dplyr::select(UUID, TPR, FPR)
cg_res <- data %>%
  filter(METHOD == "DESeq2 (CG)" & FDR == 0.05) %>%
  dplyr::select(UUID, TPR, FPR)

combined <- nocg_res %>%
  left_join(cg_res, by = "UUID") %>%
  mutate(spec_nocg = 1 - FPR.x,
         spec_cg = 1 - FPR.y) %>%
  mutate(delta_sens = TPR.y - TPR.x,
         delta_spec = spec_cg - spec_nocg)

pbottom_a <- ggplot(combined, aes(x = delta_sens)) +
  geom_histogram(bins = 30, color = "white") +
  theme_bw() +
  xlim(c(-0.5, 1)) +
  labs(x = "change in DESeq2 sensitivity\n(control genes vs. baseline)")

pbottom_b <- ggplot(combined, aes(x = delta_spec)) +
  geom_histogram(bins = 30, color = "white") +
  theme_bw() +
  xlim(c(-0.5, 1)) +
  labs(x = "change in DESeq2 specificity\n(control genes vs. baseline)")

pbottom <- plot_grid(pbottom_a, pbottom_b, ncol = 2,
                     labels = c("c", "d"),
                     scale = 0.95,
                     label_size = 18, label_y = 1.02)

p <- plot_grid(ptop, NULL, pbottom, ncol = 1,
               rel_heights = c(1, 0.05, 0.75))

ggsave(file.path(output_dir, "control_genes.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 9,
       width = 8)

# ------------------------------------------------------------------------------
#   Box plots: false positive count vs. % differential features
# ------------------------------------------------------------------------------

edgeR_only <- FALSE

# Calculate the absolute counts of false positives
if(!exists("fp")) {
  fp <- numeric(nrow(data))
  for(i in 1:nrow(data)) {
    if(i %% 10000 == 0) {
      cat(paste0("Iteration ", i, "\n"))
    }
    x <- as.numeric(str_split(data[i,]$ORACLE_BASELINE, ";")[[1]])
    y <- as.numeric(str_split(data[i,]$CALLS, ";")[[1]])
    
    x <- p.adjust(x, method = "BH")
    y <- p.adjust(y, method = "BH")
    
    x_d <- factor(x < 0.001)
    levels(x_d) <- c("not DE", "DE")
    y_d <- factor(y < 0.001)
    levels(y_d) <- c("not DE", "DE")
    
    fp[i] <- sum(y_d == "DE" & x_d == "not DE")
  }
  data$FP_n <- fp
}

if(edgeR_only) {
  temp2 <- data %>%
    filter(METHOD %in% c("edgeR", "edgeR (TMM)")) %>%
    filter(FDR == 0.05)
  temp2$METHOD <- factor(temp2$METHOD, levels = c("edgeR (TMM)", "edgeR"))
} else {
  temp2 <- data %>%
    filter(METHOD %in% method_list) %>%
    filter(FDR == 0.05)
}
temp2$P <- paste0(temp2$P, " features")
temp2 <- temp2 %>%
  mutate(PDIFF_DISCRETE = case_when(
    PERCENT_DIFF_REALIZ < 0.2 ~ "0-20%",
    PERCENT_DIFF_REALIZ < 0.4 ~ "20-40%",
    PERCENT_DIFF_REALIZ < 0.6 ~ "40-60%",
    PERCENT_DIFF_REALIZ < 0.8 ~ "60-80%",
    TRUE ~ "80-100%"
  ))
temp2$PDIFF_DISCRETE <- factor(temp2$PDIFF_DISCRETE,
                               levels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%"))

temp2$title <- paste0(temp2$METHOD, " x ", temp2$P)
p1 <- ggplot(temp2 %>% filter(P == "100 features"), aes(x = PDIFF_DISCRETE, y = FP_n, fill = METHOD)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(. ~ title, ncol = 1) + #, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = palette) +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  labs(x = "\npercent differentially abundant features",
       y = "false positive count\n") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
p2 <- ggplot(temp2 %>% filter(P == "1000 features"), aes(x = PDIFF_DISCRETE, y = FP_n, fill = METHOD)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(. ~ title, ncol = 1) + #, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = palette) +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  labs(x = "\npercent differentially abundant features",
       y = "false positive count\n") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
p3 <- ggplot(temp2 %>% filter(P == "5000 features"), aes(x = PDIFF_DISCRETE, y = FP_n, fill = METHOD)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(. ~ title, ncol = 1) + #, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = palette) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(x = "\npercent differentially abundant features",
       y = "false positive count\n") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

p1_padded <- plot_grid(p1, NULL, ncol = 1, rel_heights = c(1, ifelse(edgeR_only, 0.1, 0.05)))
p3_padded <- plot_grid(p3, NULL, ncol = 1, rel_heights = c(1, ifelse(edgeR_only, 0.1, 0.05)))
p <- plot_grid(p1_padded, p2, p3_padded, ncol = 3, rel_widths = c(1, 0.85, 0.85))

ggsave(file.path(output_dir, ifelse(edgeR_only, "F3_edgeR.svg", "F3.svg")),
       p,
       dpi = 100,
       units = "in",
       height = ifelse(edgeR_only, 4, 8),
       width = 8)

# ------------------------------------------------------------------------------
#   Box plots: false positive count vs. % differential features (REAL DATA)
# ------------------------------------------------------------------------------

dpalette <- colorRampPalette(brewer.pal(12, "Paired"))(12)

real_calls <- NULL
real_files <- list.files(path = "output/real_data_calls/no_norm", pattern = ".*\\.rds", full.names = TRUE)
for(file in real_files) {
  calls <- readRDS(file)
  file <- str_replace(file, "edgeR_TMM", "edgeR (TMM)")
  pieces <- str_split(str_split(str_split(file, "\\/")[[1]][4], "\\.")[[1]][1], "_")[[1]]
  real_calls <- rbind(real_calls,
                      data.frame(method = pieces[3],
                                 dataset = pieces[4],
                                 tpr = calls$rates$TPR,
                                 fpr = calls$rates$FPR,
                                 fp = sum(calls$rates$FP_calls),
                                 tp = sum(calls$rates$TP_calls),
                                 true_diff = sum(p.adjust(calls$all_calls$oracle_calls, method = "BH") < 0.05),
                                 n = length(calls$all_calls$oracle_calls)))
}

temp <- real_calls %>%
  mutate(true_diff_percent = true_diff/n) %>%
  arrange(true_diff_percent) %>%
  mutate(index = 1:n())

temp <- temp %>%
  mutate(datatype = case_when(
    dataset %in% c("Barlow", "VieiraSilva") ~ "16S",
    dataset %in% c("Hagai", "Monaco", "Song") ~ "bulk",
    TRUE ~ "sc"
  ))

dummy_df <- data.frame(x = c(1,2,3,4), y = c(1,2,3,4), method = c("ALDEx2", "DESeq2", "edgeR (TMM)", "scran"))
p0 <- ggplot(dummy_df, aes(x = x, y = y, shape = method)) +
  geom_point(fill = "black") +
  scale_shape_manual(values = c(21,22,24,25)) +
  theme(legend.position = "bottom")
shape_legend <- get_legend(p0)

point_sz <- 4
jitter <- 0
p <- ggplot() +
  geom_jitter(data = temp %>% filter(method == "ALDEx2"),
              mapping = aes(x = true_diff_percent*100, y = fp/n*100, fill = dataset),
              size = point_sz,
              shape = 21,
              width = jitter) +
  scale_fill_manual(values = dpalette) +
  theme(legend.position = "bottom")
legend <- get_legend(p)
p <- p +
  geom_jitter(data = temp %>% filter(method == "DESeq2"),
              mapping = aes(x = true_diff_percent*100, y = fp/n*100, fill = dataset),
              size = point_sz - 0.5,
              shape = 22,
              width = jitter) +
  theme(legend.position = "none")
p <- p +
  geom_jitter(data = temp %>% filter(method == "edgeR (TMM)"),
              mapping = aes(x = true_diff_percent*100, y = fp/n*100, fill = dataset),
              size = point_sz - 0.5,
              shape = 24,
              width = jitter) +
  theme(legend.position = "none")
p <- p +
  geom_jitter(data = temp %>% filter(method == "scran"),
              mapping = aes(x = true_diff_percent*100, y = fp/n*100, fill = dataset),
              size = point_sz,
              shape = 25,
              width = jitter) +
  theme_bw() +
  labs(x = "percent differential features",
       y = "false positives as percent total features") +
  theme(legend.position = "none") +
  xlim(c(0, 100))

p_allreal <- plot_grid(p, legend, shape_legend, ncol = 1, rel_heights = c(1, 0.2, 0.1))

ggsave(file.path(output_dir, "F8.svg"),
       p_allreal,
       dpi = 100,
       units = "in",
       height = 7,
       width = 7)

# ------------------------------------------------------------------------------
#
#   OLDER PLOTS
#
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#   Render real data outcomes on top of simulated data
# ------------------------------------------------------------------------------

DE_methods <- c("ALDEx2", "DESeq2", "scran")
datasets <- sort(c("Hagai", "Monaco", "Song", "Hashimshony", "Barlow", "Gruen",
              "Muraro", "Owens", "VieiraSilva", "Kimmerling", "Yu", "Klein"))
thresholds <- rep(1, length(datasets))
realdata <- NULL
for(i in 1:length(datasets)) {
  for(j in 1:length(DE_methods)) {
    temp <- readRDS(file.path("output", "real_data_calls", "no_norm",
                              paste0("calls_oracle_", DE_methods[j], "_", datasets[i], "_threshold", thresholds[i], ".rds")))
    realdata <- rbind(realdata,
                      data.frame(dataset = datasets[i],
                                 method = DE_methods[j],
                                 tpr = temp$rates$TPR,
                                 fpr = 1 - temp$rates$FPR))
  }
}
realdata$dataset <- factor(realdata$dataset, levels = datasets)

plots <- list()
for(i in 1:length(DE_methods)) {
  pl <- ggplot() +
    geom_point(data = data %>% filter(METHOD == DE_methods[i]),
               mapping = aes(x = FPR, y = TPR), color = "#bbbbbb", alpha = 0.5) +
    geom_point(data = realdata %>% filter(method == DE_methods[i]),
               mapping = aes(x = fpr, y = tpr, fill = dataset),
               shape = 21,
               size = 3,
               color = "black") +
    scale_fill_brewer(palette = "Paired") +
    xlim(c(-0.025,1.025)) +
    ylim(c(-0.025,1.025)) +
    labs(x = paste0(DE_methods[i], " specificity (1 - FPR)"),
         y = paste0(DE_methods[i], " sensitivity (TPR)")) +
    theme_bw() +
    theme(axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)))
  # pl <- ggMarginal(pl, margins = "both", type = "histogram", bins = 30,
  #                  fill = "#888888", color = "#333333", size = 10)
  if(i < 3) {
    pl <- pl +
      theme(legend.position = "none")
  }
  plots[[length(plots)+1]] <- pl
}

pl <- plot_grid(plotlist = plots, ncol = 3, rel_widths = c(1, 1, 1.3))
ggsave(file.path("output", "images", "real_results_over_simulated.png"),
       pl,
       dpi = 100,
       units = "in",
       height = 4,
       width = 14)

# ------------------------------------------------------------------------------
#   Sensitivity x specificity plot labeled for majority differential features
#   PLUS high or very high fold change
# ------------------------------------------------------------------------------

# p1_100 <- NULL; p1_1000 <- NULL; p1_5000 <- NULL
# p2_100 <- NULL; p2_1000 <- NULL; p2_5000 <- NULL
# p3_100 <- NULL; p3_1000 <- NULL; p3_5000 <- NULL
# legend <- NULL
# for(i in 1:length(method_list)) {
#   method <- method_list[i]
#   for(j in c(100, 1000, 5000)) {
#     # pl <- plot_ROC_flag(data, method, j, data$PERCENT_DIFF_REALIZ < 0.5 & data$FC_plot %in% levels(data$FC_plot)[3:4])
#     plot_pieces <- plot_ROC_fc(data %>% filter(PERCENT_DIFF_REALIZ < 0.5 & FC_plot %in% levels(data$FC_plot)[3:4]), method, j)
#     pl <- plot_pieces$pl
#     if(is.null(legend)) {
#       legend <- plot_pieces$legend
#     }
#     if(i == 1 & j == 100) p1_100 <<- pl
#     if(i == 1 & j == 1000) p1_1000 <<- pl
#     if(i == 1 & j == 5000) p1_5000 <<- pl
#     if(i == 2 & j == 100) p2_100 <<- pl
#     if(i == 2 & j == 1000) p2_1000 <<- pl
#     if(i == 2 & j == 5000) p2_5000 <<- pl
#     if(i == 3 & j == 100) p3_100 <<- pl
#     if(i == 3 & j == 1000) p3_1000 <<- pl
#     if(i == 3 & j == 5000) p3_5000 <<- pl
#   }
# }
# prow1 <- plot_grid(p1_100, p2_100, p3_100, ncol = 3)
# prow2 <- plot_grid(p1_1000, p2_1000, p3_1000, ncol = 3)
# prow3 <- plot_grid(p1_5000, p2_5000, p3_5000, ncol = 3)
# # pl <- plot_grid(prow1, prow2, prow3, nrow = 3, labels = c("a", "b", "c"), label_size = 18, label_y = 1.02)
# pgrid <- plot_grid(prow1, prow2, prow3, nrow = 3, labels = c("a", "b", "c"), label_size = 18, label_y = 1.02)
# pl <- plot_grid(pgrid, legend, ncol = 1, rel_heights = c(1, .1))
# # ggsave(file.path("output", "images", "ROC_by_subset1.png"),
# ggsave(file.path("output", "images", "ROC_by_subset2.png"),
#        plot = pl,
#        units = "in",
#        height = 10,
#        width = 10,
#        bg = "white")
# # show(pl)
# 
# temp <- data %>%
#   filter(METHOD == "scran" & FC_plot %in% levels(data$FC_plot)[4] & P == 5000)
# median(temp$FPR)

# ------------------------------------------------------------------------------
#   Sensitivity x specificity plot using discrepancies against calls made by a
#   NB GLM ("oracle method") on absolute abundances
# ------------------------------------------------------------------------------

# data <- pull_data(qq, use_baseline = "oracle")
# data$FPR <- 1 - data$FPR
# 
# p1_100 <- NULL; p1_1000 <- NULL; p1_5000 <- NULL
# p2_100 <- NULL; p2_1000 <- NULL; p2_5000 <- NULL
# p3_100 <- NULL; p3_1000 <- NULL; p3_5000 <- NULL
# legend <- NULL
# for(i in 1:length(method_list)) {
#   method <- method_list[i]
#   for(j in c(100, 1000, 5000)) {
#     # plot_pieces <- plot_ROC_fc(data, method, j)
#     plot_pieces <- plot_ROC_percentDE(data, method, j)
#     pl <- plot_pieces$pl
#     if(is.null(legend)) {
#       legend <- plot_pieces$legend
#     }
#     if(i == 1 & j == 100) p1_100 <<- pl
#     if(i == 1 & j == 1000) p1_1000 <<- pl
#     if(i == 1 & j == 5000) p1_5000 <<- pl
#     if(i == 2 & j == 100) p2_100 <<- pl
#     if(i == 2 & j == 1000) p2_1000 <<- pl
#     if(i == 2 & j == 5000) p2_5000 <<- pl
#     if(i == 3 & j == 100) p3_100 <<- pl
#     if(i == 3 & j == 1000) p3_1000 <<- pl
#     if(i == 3 & j == 5000) p3_5000 <<- pl
#   }
# }
# prow1 <- plot_grid(p1_100, p2_100, p3_100, ncol = 3)
# prow2 <- plot_grid(p1_1000, p2_1000, p3_1000, ncol = 3)
# prow3 <- plot_grid(p1_5000, p2_5000, p3_5000, ncol = 3)
# pgrid <- plot_grid(prow1, prow2, prow3, nrow = 3, labels = c("a", "b", "c"), label_size = 18, label_y = 1.02)
# pl <- plot_grid(pgrid, legend, ncol = 1, rel_heights = c(1, .1))
# ggsave(file.path("output", "images", "ROC_by_FC_NBGLM.png"),
#        plot = pl,
#        units = "in",
#        height = 10.5,
#        width = 10,
#        bg = "white")
# show(pl)

# ------------------------------------------------------------------------
#   What proportion of calls exceed an FDR of 5% here?
# ------------------------------------------------------------------------

# res <- res %>%
#   filter(METHOD == 'DESeq2') %>%
#   mutate(FDR_big = ifelse(FPR > 0.05, "yes", "no")) %>%
#   select(!c("ORACLE_BASELINE", "SELF_BASELINE", "CALLS"))
# 
# sum(res$FDR_big == "yes")/nrow(res)

# Filter out what DESeq2 calls > 85% DE
# res$DESeq2_DE <- sapply(res$SELF_BASELINE, function(x) sum(string_to_calls(x)))
# res$DESeq2_DE <- res$DESeq2_DE / p
# res <- res %>%
#   filter(DESeq2_DE <= 0.85)

# res <- res %>%
#   filter(PERCENT_DIFF_REALIZ < 0.85)

# ------------------------------------------------------------------------
#   Histogram with real data overlaid
# ------------------------------------------------------------------------

# real_data <- readRDS("bad_good_results.rds")
# real_data <- real_data %>%
#   filter(dtype != "simulated") %>%
#   select(dtype, fpr)
# 
# real_df <- data.frame(x = real_data$fpr,
#                       y = rep(c(25,45), 4),
#                       type = real_data$dtype)
# ggplot() +
#   geom_histogram(data = res[res$METHOD == "DESeq2",],
#                  mapping = aes(x = FPR),
#                  color = "white") +
#   geom_point(data = real_df,
#              mapping = aes(x = x, y = y, fill = type),
#              size = 4,
#              shape = 21,
#              stroke = 1.2) +
#   scale_fill_manual(values = generate_highcontrast_palette(8)) +
#   labs(fill = "Data set") +
#   theme_bw()

# Remove some of the simulations with huge median abundances
# res <- res %>%
#   filter(MED_ABS_TOTAL < 2e06)

# Remove simulations where > 50% of features are differential
# res <- res %>%
#   filter(PERCENT_DIFF < 0.5)