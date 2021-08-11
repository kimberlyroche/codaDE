library(tidyverse)
library(codaDE)
library(rulesoflife) # for generate_highcontrast_palette()

palette <- generate_highcontrast_palette(5001)

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
    scale_fill_manual(values = palette)
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
  feat_idx <- data.frame(x = 1:nrow(relab), mu = rowMeans(relab)) %>%
    arrange(desc(mu)) %>%
    top_n(k) %>%
    arrange(x) %>%
    pull(x)
  relab1 <- relab[feat_idx,]
  relab2 <- colSums(relab[-feat_idx,])
  relab <- rbind(relab1, relab2)
  
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

  p <- ggplot(data_long, aes(fill = feature, y = abundance, x = sample)) + 
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = palette) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
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

get_feature_count <- function(counts) {
  nrow(counts[rowMeans(counts) >= 1,])
}

# ------------------------------------------------------------------------------
#   Vieira-Silva et al. data
# ------------------------------------------------------------------------------

use_groups <- c("mHC", "CD")
vs <- parse_VieiraSilva(absolute = TRUE)
cat(paste0("Feature count: ", get_feature_count(vs$counts), "\n"))
plot_PCA(vs$counts, vs$groups, use_groups = use_groups, save_name = "vieirasilva_01.png")
calc_prop_DA(vs$counts, vs$groups, use_groups = use_groups)
calc_FC(vs$counts, vs$groups, use_groups = use_groups)
plot_totals(vs$counts, vs$groups, use_groups = use_groups, save_name = "vieirasilva_02.png")
# Subset to CD vs. mHC samples
idx_subset <- c(sample(which(vs$groups == "mHC"), size = 20),
                sample(which(vs$groups == "CD"), size = 20))
plot_relab(vs$counts[,idx_subset], vs$groups[idx_subset], use_props = FALSE, save_name = "vieirasilva_03a.png")
plot_relab(vs$counts[,idx_subset], vs$groups[idx_subset], use_props = TRUE, save_name = "vieirasilva_03b.png")

# ------------------------------------------------------------------------------
#   Barlow et al. data
# ------------------------------------------------------------------------------

barlow <- parse_Barlow(absolute = TRUE)
cat(paste0("Feature count: ", get_feature_count(barlow$counts), "\n"))
plot_PCA(barlow$counts, barlow$groups, save_name = "barlow_01.png")
calc_prop_DA(barlow$counts, barlow$groups)
calc_FC(barlow$counts, barlow$groups)
plot_totals(barlow$counts, barlow$groups, save_name = "barlow_02.png")
plot_relab(barlow$counts, barlow$groups, use_props = FALSE, save_name = "barlow_03a.png")
plot_relab(barlow$counts, barlow$groups, use_props = TRUE, save_name = "barlow_03b.png")

# ------------------------------------------------------------------------------
#   Song et al. data
# ------------------------------------------------------------------------------

song <- parse_Song(absolute = TRUE)
cat(paste0("Feature count: ", get_feature_count(song$counts), "\n"))
plot_PCA(song$counts, song$groups, save_name = "song_01.png")
calc_prop_DA(song$counts, song$groups)
calc_FC(song$counts, song$groups)
plot_totals(song$counts, song$groups, save_name = "song_02.png")
plot_relab(song$counts, song$groups, k = 200, use_props = FALSE, save_name = "song_03a.png")
plot_relab(song$counts, song$groups, k = 200, use_props = TRUE, save_name = "song_03b.png")

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

use_groups <- c("CD4_naive", "PBMC")
monaco <- parse_Monaco(absolute = TRUE)
cat(paste0("Feature count: ", get_feature_count(monaco$counts), "\n"))
plot_PCA(monaco$counts, monaco$groups, use_groups = use_groups,
         use_labels = TRUE, save_name = "monaco_01.png", height = 6, width = 9)
calc_prop_DA(monaco$counts, monaco$groups, use_groups = use_groups)
calc_FC(monaco$counts, monaco$groups, use_groups = use_groups)
plot_totals(monaco$counts, monaco$groups, use_groups = use_groups,
            save_name = "monaco_02.png", height = 6, width = 10)
plot_relab(monaco$counts, monaco$group, use_groups = use_groups, k = 5000,
           use_props = FALSE, save_name = "monaco_03a.png", height = 6, width = 12)
plot_relab(monaco$counts, monaco$group, use_groups = use_groups, k = 5000,
           use_props = TRUE, save_name = "monaco_03b.png", height = 6, width = 12)

# ------------------------------------------------------------------------------
#   Muraro et al. data
# ------------------------------------------------------------------------------

muraro <- parse_Muraro(absolute = TRUE)
cat(paste0("Feature count: ", get_feature_count(muraro$counts), "\n"))
idx_subset <- c(sample(which(muraro$groups == "alpha"), size = 50),
                sample(which(muraro$groups == "beta"), size = 50))
muraro$counts <- muraro$counts[,idx_subset]
muraro$groups <- muraro$groups[idx_subset]
plot_PCA(muraro$counts, muraro$groups, save_name = "muraro_01.png")
calc_prop_DA(muraro$counts, muraro$groups)
calc_FC(muraro$counts, muraro$groups)
# Subset the groups for the following
plot_totals(muraro$counts, muraro$groups, save_name = "muraro_02.png")
plot_relab(muraro$counts, muraro$groups, k = 5000, use_props = FALSE, save_name = "muraro_03a.png")
plot_relab(muraro$counts, muraro$groups, k = 5000, use_props = TRUE, save_name = "muraro_03b.png")

# ------------------------------------------------------------------------------
#   Hashimshony et al. data
# ------------------------------------------------------------------------------

use_groups = c("0", "0.5")
hashim <- parse_Hashimshony(absolute = TRUE)
cat(paste0("Feature count: ", get_feature_count(hashim$counts), "\n"))
plot_PCA(hashim$counts, hashim$groups, save_name = "hashimshony_01.png")
calc_prop_DA(hashim$counts, hashim$groups, use_groups = use_groups)
calc_FC(hashim$counts, hashim$groups, use_groups = use_groups)
plot_totals(hashim$counts, hashim$groups, save_name = "hashimshony_02.png")
plot_relab(hashim$counts, hashim$groups, use_props = FALSE, use_groups = use_groups,
           k = 5000, save_name = "hashimshony_03a.png")
plot_relab(hashim$counts, hashim$groups, use_props = TRUE, use_groups = use_groups,
           k = 5000, save_name = "hashimshony_03b.png")

# ------------------------------------------------------------------------------
#   Kimmerling et al. data
# ------------------------------------------------------------------------------

kimmer <- parse_Kimmerling(absolute = TRUE)
idx_subset <- c(sample(which(kimmer$groups == "low_mass"), size = 20),
                sample(which(kimmer$groups == "high_mass"), size = 20))
kimmer$counts <- kimmer$counts[,idx_subset]
kimmer$groups <- kimmer$groups[idx_subset]
cat(paste0("Feature count: ", get_feature_count(kimmer$counts), "\n"))
plot_PCA(kimmer$counts, kimmer$groups, save_name = "kimmerling_01.png")
calc_prop_DA(kimmer$counts, kimmer$groups)
calc_FC(kimmer$counts, kimmer$groups)
# Subset the groups for the following
plot_totals(kimmer$counts, kimmer$groups, save_name = "kimmerling_02.png")
plot_relab(kimmer$counts, kimmer$group, use_props = FALSE, k = 5000,
           save_name = "kimmerling_03a.png")
plot_relab(kimmer$counts, kimmer$group, use_props = TRUE, k = 5000,
           save_name = "kimmerling_03b.png")
















