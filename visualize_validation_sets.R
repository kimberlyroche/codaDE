library(tidyverse)
library(codaDE)
library(rulesoflife) # for generate_highcontrast_palette()

img_dir <- file.path("output", "images", "validation_datasets")

# `counts` should be features x samples
plot_PCA <- function(counts, groups, use_labels = FALSE, save_name = NULL, height = 4, width = 6) {
  if(nrow(counts) > 1000) {
    counts <- counts[sample(1:nrow(counts), size = 1000),]
  }
  coords <- cmdscale(dist(t(log(counts + 0.5))))
  palette <- generate_highcontrast_palette(length(unique(groups)))
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
plot_totals <- function(counts, groups, save_name = NULL, height = 4, width = 6) {
  plot_df <- data.frame(totals = colSums(counts),
                        label = groups)
  plot_df <- plot_df %>%
    arrange(label)
  plot_df$x <- 1:nrow(plot_df)
  palette <- generate_highcontrast_palette(length(unique(groups)))
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
plot_relab <- function(counts, groups, k = 20, save_name = NULL, height = 4, width = 6) {
  # Subset to top k=20 features
  relab <- counts
  for(j in 1:ncol(relab)) {
    relab[,j] <- relab[,j] / sum(relab[,j])
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
  
  palette <- generate_highcontrast_palette(length(unique(data_long$feature)))
  
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

# ------------------------------------------------------------------------------
#   Barlow et al. data
# ------------------------------------------------------------------------------

barlow <- parse_Barlow(absolute = TRUE)
nrow(barlow$counts[rowMeans(barlow$counts) >= 1,])
plot_PCA(barlow$counts, barlow$groups, save_name = "barlow_01.png")
calc_prop_DA(barlow$counts, barlow$groups)
calc_FC(barlow$counts, barlow$groups)
plot_totals(barlow$counts, barlow$groups, save_name = "barlow_02.png")
plot_relab(barlow$counts, barlow$groups, save_name = "barlow_03.png")

# ------------------------------------------------------------------------------
#   Morton et al. data
#
#   Excluding this data set. It's hugely noisy with essentially no detectable
#   differential abundance.
# ------------------------------------------------------------------------------

# morton <- parse_Morton(absolute = TRUE)
# nrow(morton$counts[rowMeans(morton$counts) >= 1,])
# plot_PCA(morton$counts, morton$groups, save_name = "morton_01.png")
# calc_prop_DA(morton$counts, morton$groups)
# calc_FC(morton$counts, morton$groups)
# plot_totals(morton$counts, morton$groups, save_name = "morton_02.png")
# plot_relab(morton$counts, morton$groups, save_name = "morton_03.png")

# ------------------------------------------------------------------------------
#   Song et al. data
# ------------------------------------------------------------------------------

song <- parse_Song(absolute = TRUE)
nrow(song$counts[rowMeans(song$counts) >= 1,])
plot_PCA(song$counts, song$groups, save_name = "song_01.png")
calc_prop_DA(song$counts, song$groups)
calc_FC(song$counts, song$groups)
plot_totals(song$counts, song$groups, save_name = "song_02.png")
plot_relab(song$counts, song$groups, k = 200, save_name = "song_03.png")

# ------------------------------------------------------------------------------
#   Muraro et al. data
# ------------------------------------------------------------------------------

muraro <- parse_Muraro(absolute = TRUE)
nrow(muraro$counts[rowMeans(muraro$counts) >= 1,])
plot_PCA(muraro$counts, muraro$groups, save_name = "muraro_01.png")
calc_prop_DA(muraro$counts, muraro$groups) # debug SF problem
calc_FC(muraro$counts, muraro$groups)
# Subset the groups for the following
idx_subset <- c(sample(which(muraro$groups == "alpha"), size = 20),
                sample(which(muraro$groups == "beta"), size = 20))
plot_totals(muraro$counts[,idx_subset], muraro$groups[idx_subset], save_name = "muraro_02.png")
plot_relab(muraro$counts[,idx_subset], muraro$groups[idx_subset], k = 5000, save_name = "muraro_03.png")

# ------------------------------------------------------------------------------
#   Monaco et al. data
# ------------------------------------------------------------------------------

monaco <- parse_Monaco(absolute = TRUE)
nrow(monaco$counts[rowMeans(monaco$counts) >= 1,])
plot_PCA(monaco$counts, monaco$groups, use_labels = TRUE, save_name = "monaco_01.png", height = 6, width = 9)
calc_prop_DA(monaco$counts, monaco$groups, use_groups = c("Plasmablasts", "Neutrophils"))
# calc_prop_DA(monaco$counts, monaco$groups, use_groups = c("PBMC", "CD4_naive"))
calc_FC(monaco$counts, monaco$groups, use_groups = c("Plasmablasts", "Neutrophils"))
plot_totals(monaco$counts, monaco$groups, save_name = "monaco_02.png", height = 6, width = 10)
plot_relab(monaco$counts, monaco$group, k = 5000, save_name = "monaco_03.png", height = 6, width = 12)

# ------------------------------------------------------------------------------
#   Hashimshony et al. data
# ------------------------------------------------------------------------------

hashim <- parse_Hashimshony(absolute = TRUE, use_spike_ins = FALSE)
nrow(hashim$counts[rowMeans(hashim$counts) >= 1,])
plot_PCA(hashim$counts, hashim$groups, save_name = "hashimshony_01.png")
calc_prop_DA(hashim$counts, hashim$groups, use_groups = c("0", "1"))
# calc_prop_DA(monaco$counts, monaco$groups, use_groups = c("PBMC", "CD4_naive"))
calc_FC(hashim$counts, hashim$groups, use_groups = c("0", "1"))
plot_totals(hashim$counts, hashim$groups, save_name = "hashimshony_02.png")
plot_relab(hashim$counts, hashim$group, k = 5000, save_name = "hashimshony_03.png")

# ------------------------------------------------------------------------------
#   Kimmerling et al. data
# ------------------------------------------------------------------------------

kimmer <- parse_Kimmerling(absolute = TRUE)
nrow(kimmer$counts[rowMeans(kimmer$counts) >= 1,])
plot_PCA(kimmer$counts, kimmer$groups, save_name = "kimmerling_01.png")
calc_prop_DA(kimmer$counts, kimmer$groups)
calc_FC(kimmer$counts, kimmer$groups)
# Subset the groups for the following
idx_subset <- c(sample(which(kimmer$groups == "low_mass"), size = 20),
                sample(which(kimmer$groups == "high_mass"), size = 20))
plot_totals(kimmer$counts[,idx_subset], kimmer$groups[idx_subset], save_name = "kimmerling_02.png")
plot_relab(kimmer$counts[,idx_subset], kimmer$group[idx_subset], k = 5000, save_name = "kimmerling_03.png")

# ------------------------------------------------------------------------------
#   Vieira-Silva et al. data
# ------------------------------------------------------------------------------

vs <- parse_VieiraSilva(absolute = TRUE)
nrow(vs$counts[rowMeans(vs$counts) >= 1,])
plot_PCA(vs$counts, vs$groups, save_name = "vieirasilva_01.png")
calc_prop_DA(vs$counts, vs$groups, use_groups = c("mHC", "CD"))
calc_FC(vs$counts, vs$groups, use_groups = c("mHC", "CD"))
plot_totals(vs$counts, vs$groups, save_name = "vieirasilva_02.png")
# Subset to CD vs. mHC samples
idx_subset <- c(sample(which(vs$groups == "mHC"), size = 20),
                sample(which(vs$groups == "CD"), size = 20))
plot_relab(vs$counts[,idx_subset], vs$groups[idx_subset], save_name = "vieirasilva_03.png")
