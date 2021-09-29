library(tidyverse)
library(codaDE)
library(rulesoflife) # for generate_highcontrast_palette()

source("ggplot_fix.R")

palette <- generate_highcontrast_palette(6000)

img_dir <- file.path("output", "images", "validation_datasets")

# `counts` should be features x samples
plot_PCA <- function(counts, groups, use_labels = FALSE, use_groups = NULL,
                     save_name = NULL, height = 4, width = 6) {
  if(nrow(counts) > 1000) {
    counts <- counts[sample(1:nrow(counts), size = 1000),]
  }
  if(!is.null(use_groups)) {
    counts <- counts[,groups %in% use_groups]
    groups <- as.character(groups[groups %in% use_groups])
  }
  coords <- cmdscale(dist(t(log(counts + 0.5))))
  p <- ggplot(data.frame(x = coords[,1], y = coords[,2], label = groups),
              aes(x = x, y = y, fill = label, label = label)) +
    geom_point(size = 3, shape = 21) +
    theme_bw() +
    labs(x = "PC 1", y = "PC 2", fill = "Condition") +
    scale_fill_manual(values = c("#e6780b", "#dddddd"))
  if(use_labels) {
    p <- p +
      geom_text(check_overlap = TRUE)
  }
  show(p)
  if(!is.null(save_name)) {
    ggsave(file.path(img_dir, save_name), p,
           units = "in", height = height, width = width, dpi = 100)
  }
}

# `counts` should be features x samples
plot_totals <- function(counts, groups, use_groups = NULL, save_name = NULL,
                        height = 4, width = 6) {
  if(!is.null(use_groups)) {
    counts <- counts[,groups %in% use_groups]
    groups <- as.character(groups[groups %in% use_groups])
  }
  plot_df <- data.frame(totals = colSums(counts),
                        label = groups)
  plot_df <- plot_df %>%
    arrange(label)
  plot_df$x <- 1:nrow(plot_df)
  p <- ggplot(plot_df, aes(x = x, y = totals, fill = label)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    labs(x = "sample index", y = "total abundance", fill = "Condition") +
    scale_fill_manual(values = palette)
  show(p)
  if(!is.null(save_name)) {
    ggsave(file.path(img_dir, save_name), p,
           units = "in", height = height, width = width, dpi = 100)
  }
}

# `counts` should be features x samples
plot_relab <- function(counts, groups, k = 20, use_groups = NULL, use_props = TRUE,
                       save_name = NULL, height = 4, width = 6) {
  if(!is.null(use_groups)) {
    counts <- counts[,groups %in% use_groups]
    groups <- as.character(groups[groups %in% use_groups])
  }
  # Subset to top k=20 features
  relab <- counts
  if(use_props) {
    for(j in 1:ncol(relab)) {
      relab[,j] <- relab[,j] / sum(relab[,j])
    }
  }
  
  # START -- Comment out this section to force sample palettes
  feat_idx <- suppressMessages(data.frame(x = 1:nrow(relab), mu = rowMeans(relab)) %>%
                                 arrange(desc(mu)) %>%
                                 top_n(k) %>%
                                 arrange(x) %>%
                                 pull(x))
  relab1 <- relab[feat_idx,]
  relab2 <- colSums(relab[-feat_idx,])
  relab <- rbind(relab1, relab2)
  # END -- Comment out
  
  # `relab` is oriented as FEATURES x SAMPLES
  n_samples <- ncol(relab)
  relab <- cbind(1:nrow(relab), relab)
  colnames(relab) <- c("feature", paste0(groups, " (sample ", 1:length(groups), ")"))
  data_long <- pivot_longer(as.data.frame(relab),
                            !feature,
                            names_to = "sample",
                            values_to = "abundance")

  data_long$feature <- factor(data_long$feature)
  data_long$sample <- factor(data_long$sample)

  tick_labels <- as.character(groups)
  names(tick_labels) <- colnames(relab)[2:ncol(relab)]
  p <- ggplot(data_long, aes(fill = feature, y = abundance, x = sample)) + 
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = palette) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_x_discrete(labels = tick_labels)
  show(p)
  if(!is.null(save_name)) {
    ggsave(file.path(img_dir, save_name), p,
           units = "in", height = height, width = width, dpi = 100)
  }
}

# `counts` should be passed as features x samples
calc_prop_DA <- function(counts, groups, use_groups = NULL) {
  # Filter lowly abundant stuff
  counts <- counts[rowMeans(counts) >= 1,]
  if(nrow(counts) > 1000) {
    counts <- counts[sample(1:nrow(counts), size = 1000),]
  }
  if(!is.null(use_groups)) {
    counts <- counts[,groups %in% use_groups]
    groups <- as.character(groups[groups %in% use_groups])
  }
  # Percent differentially expressed features
  rounded_counts <- t(round(counts))
  res <- call_DA_DESeq2(rounded_counts, groups)
  x <- length(res$pval)
  y <- sum(p.adjust(res$pval, method = "BH") < 0.05)
  cat(paste0("Differential: ", round(y, 3), " / ", x, " (", round(y/x, 3)*100, "%)\n"))
}

# `counts` should be passed as features x samples
calc_FC <- function(counts, groups, use_groups = NULL) {
  if(!is.null(use_groups)) {
    counts <- counts[,groups %in% use_groups]
    groups <- groups[groups %in% use_groups]
  }
  counts_A <- counts[,groups == unique(groups)[1]]
  counts_B <- counts[,groups == unique(groups)[2]]
  x <- mean(colSums(counts_A))
  y <- mean(colSums(counts_B))
  ratio <- max(x, y) / min(x, y)
  cat(paste0("Absolute FC: ", round(ratio, 3), "\n"))
}

get_feature_count <- function(counts, threshold = 1) {
  nrow(counts[rowMeans(counts) >= threshold,])
}

get_zeros <- function(counts, threshold = 1) {
  retain_features <- colMeans(counts) >= threshold & colMeans(counts) >= threshold
  counts <- counts[,retain_features]
  sum(counts == 0)/(nrow(counts)*ncol(counts))
}

summarize_ds <- function(parse_fn, dataset_name = NULL) {
  data <- do.call(parse_fn, list(absolute = TRUE))
  counts <- data$counts
  groups <- factor(data$groups)
  
  cat(paste0("Feature count: ", get_feature_count(counts), "\n"))
  calc_FC(counts, groups)
  calc_prop_DA(counts, groups)
  cat(paste0("Zeros: ", round(get_zeros(counts, threshold = 5)*100), "%\n\n"))
  
  # plot_PCA(counts,
  #          groups,
  #          save_name = ifelse(!is.null(dataset_name), paste0(dataset_name, "_01.png"), NULL))
  g1_idx <- which(groups == unique(groups)[1])
  g2_idx <- which(groups == unique(groups)[2])
  if(length(g1_idx) > 20) {
    g1_idx <- sample(g1_idx, size = 20)
  }
  if(length(g2_idx) > 20) {
    g2_idx <- sample(g2_idx, size = 20)
  }
  idx_subset <- c(g1_idx, g2_idx)
  use_k <- min(round(nrow(counts)*0.8), 5000)
  if(is.null(dataset_name)) {
    save_name <- NULL
  } else {
    save_name <- paste0(dataset_name, "_01_abs.png")
  }
  plot_relab(counts[,idx_subset],
             groups[idx_subset],
             k = use_k,
             use_props = FALSE,
             save_name = save_name)
  
  data <- do.call(parse_fn, list(absolute = FALSE))
  counts <- data$counts
  groups <- factor(data$groups)

  if(is.null(dataset_name)) {
    save_name <- NULL
  } else {
    save_name <- paste0(dataset_name, "_02_rel.png")
  }
  plot_relab(counts[,idx_subset],
             groups[idx_subset],
             k = use_k,
             use_props = FALSE,
             save_name = save_name)
}

# ------------------------------------------------------------------------------
#   Vieira-Silva et al. data
# ------------------------------------------------------------------------------

summarize_ds("parse_VieiraSilva", dataset_name = "VieiraSilva")

# ------------------------------------------------------------------------------
#   Barlow et al. data
# ------------------------------------------------------------------------------

summarize_ds("parse_Barlow", dataset_name = "Barlow")

# ------------------------------------------------------------------------------
#   Song et al. data
# ------------------------------------------------------------------------------

summarize_ds("parse_Song", dataset_name = "Song")

# # PCA of features in FPR RF training set vs. this data set
# # First pull "features_df" from `validate.R` for Song et al. x MAST
# fit_obj <- readRDS(file.path("output", "predictive_fits", "all", "all_FPR_self.rds"))
# features2 <- fit_obj$train_features %>% filter(P == 1000 & METHOD == "MAST")
# subset_n <- 100
# features2 <- features2[sample(1:nrow(features2), size = subset_n),]
# features <- rbind(features_df, features2)
# features$METHOD <- as.numeric(features$METHOD)
# features <- apply(features, 2, scale)
# d <- dist(features)
# coords <- cmdscale(d)
# 
# # [1] "TOTALS_C_FC"    "TOTALS_C_D"     "TOTALS_C_MAX_D" "TOTALS_C_MED_D"
# # [5] "TOTALS_C_SD_D"  "CORR_RA_MED"    "CORR_RA_SD"     "CORR_RA_SKEW"  
# # [9] "CORR_LOG_MED"   "CORR_LOG_SD"    "CORR_LOG_SKEW"  "CORR_CLR_MED"  
# # [13] "CORR_CLR_SD"    "CORR_CLR_SKEW"  "COMP_C_P0_A"    "COMP_C_P0_B"   
# # [17] "COMP_C_P1_A"    "COMP_C_P1_B"    "COMP_C_P5_A"    "COMP_C_P5_B"   
# # [21] "COMP_RA_P01_A"  "COMP_RA_P01_B"  "COMP_RA_P1_A"   "COMP_RA_P1_B"  
# # [25] "COMP_RA_P5_A"   "COMP_RA_P5_B"   "COMP_RA_MAX_A"  "COMP_RA_MED_A" 
# # [29] "COMP_RA_SD_A"   "COMP_RA_SKEW_A" "COMP_RA_MAX_B"  "COMP_RA_MED_B" 
# # [33] "COMP_RA_SD_B"   "COMP_RA_SKEW_B" "COMP_C_ENT_A"   "COMP_C_ENT_B"  
# # [37] "FW_RA_MAX_D"    "FW_RA_MED_D"    "FW_RA_SD_D"     "FW_RA_PPOS_D"  
# # [41] "FW_RA_PNEG_D"   "FW_RA_PFC05_D"  "FW_RA_PFC2_D"   "FW_LOG_MAX_D"  
# # [45] "FW_LOG_MED_D"   "FW_LOG_SD_D"    "FW_LOG_PPOS_D"  "FW_LOG_PNEG_D" 
# # [49] "FW_LOG_PFC05_D" "FW_LOG_PFC1_D"  "FW_LOG_PFC2_D"  "FW_CLR_MAX_D"  
# # [53] "FW_CLR_PPOS_D"  "FW_CLR_PFC05_D" "FW_CLR_PFC1_D"  "FW_CLR_PFC2_D" 
# # [57] "P"              "METHOD"  
# 
# label_by <- "TOTALS_C_FC"
# ggplot(data.frame(x = coords[,1],
#                   y = coords[,2],
#                   identity = c("Song", rep("training", subset_n)),
#                   label = features[,colnames(features) == label_by]),
#        aes(x = x, y = y, fill = label, color = identity)) +
#   geom_point(size = 3, shape = 21, stroke = 2) +
#   scale_color_manual(values = c("red", "black")) +
#   theme_bw() +
#   labs(x = "PC 1", y = "PC 2", fill = label_by)

# ------------------------------------------------------------------------------
#   Monaco et al. data
# ------------------------------------------------------------------------------

summarize_ds("parse_Monaco", dataset_name = "Monaco")

# ------------------------------------------------------------------------------
#   Hagai et al. data
# ------------------------------------------------------------------------------

summarize_ds("parse_Hagai", dataset_name = "Hagai")

# ------------------------------------------------------------------------------
#   Owens et al. data
# ------------------------------------------------------------------------------

summarize_ds("parse_Owens", dataset_name = "Owens")

# ------------------------------------------------------------------------------
#   Klein et al. data
# ------------------------------------------------------------------------------

summarize_ds("parse_Klein", dataset_name = "Klein")

# ------------------------------------------------------------------------------
#   Yu et al. data
# ------------------------------------------------------------------------------

summarize_ds("parse_Yu", dataset_name = "Yu")

# ------------------------------------------------------------------------------
#   Unused data sets
# ------------------------------------------------------------------------------

summarize_ds("parse_Morton", dataset_name = "Morton")
summarize_ds("parse_Athanasiadou", dataset_name = "Athanasiadou")
summarize_ds("parse_Muraro", dataset_name = "Muraro")
summarize_ds("parse_Hashimshony", dataset_name = "Hashimshony")
summarize_ds("parse_Kimmerling", dataset_name = "Kimmerling")
summarize_ds("parse_ESCA", dataset_name = "ESCA")
summarize_ds("parse_Gruen", dataset_name = "Gruen")
summarize_ds("parse_Ferreira", dataset_name = "Ferreira")

























