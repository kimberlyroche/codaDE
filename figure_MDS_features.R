source("path_fix.R")

library(tidyverse)
library(codaDE)
library(ggsci)
library(ggbiplot)
#library("devtools")
#install_github("vqv/ggbiplot")

start <- Sys.time()

datasets <- c("VieiraSilva", "Barlow", "Song", "Monaco", "Hagai", "Owens", "Klein", "Yu")
thresholds <- c(1, 1, 1, 2, 1, 1, 1, 1)

save_fn <- file.path("output", "MDS_features.rds")
if(file.exists(save_fn)) {
  temp_obj <- readRDS(save_fn)
  real_features <- temp_obj$real_features
  simulated_features <- temp_obj$simulated_features
} else {
  real_features <- NULL
  for(i in 1:length(datasets)) {
    dataset_name <- datasets[i]
    threshold <- thresholds[i]

    cat(paste0("Parsing dataset: ", dataset_name, "\n"))
  
    abs_data <- do.call(paste0("parse_", dataset_name), list(absolute = TRUE))
    rel_data <- do.call(paste0("parse_", dataset_name), list(absolute = FALSE))
  
    # Reorient as (samples x features)
    ref_data <- t(abs_data$counts)
    data <- t(rel_data$counts)
    groups <- abs_data$groups
    groups <- factor(groups)
  
    # Subsample if tons of cells/samples
    downsample_limit <- 100 # was 100
    set.seed(100)
    pairs <- table(groups)
    if(pairs[1] > downsample_limit) {
      A_sample <- sample(which(groups == names(pairs)[1]), size = downsample_limit)
    } else {
      A_sample <- which(groups == names(pairs)[1])
    }
    if(pairs[2] > downsample_limit) {
      B_sample <- sample(which(groups == names(pairs)[2]), size = downsample_limit)
    } else {
      B_sample <- which(groups == names(pairs)[2])
    }
    ref_data <- ref_data[c(A_sample, B_sample),]
    data <- data[c(A_sample, B_sample),]
    groups <- groups[c(A_sample, B_sample)]

    # Convert to integers, just for DESeq2
    ref_data <- apply(ref_data, c(1,2), as.integer)
    data <- apply(data, c(1,2), as.integer)

    # Look at sparsity; filter out too-low features
    retain_features <- colMeans(ref_data) >= threshold & colMeans(data) >= threshold
    ref_data <- ref_data[,retain_features]
    data <- data[,retain_features]
  
    # Get info we need from data to make a prediction
    if(dataset_name == "VieiraSilva") {
      counts_A <- data[groups == "mHC",]
      counts_B <- data[groups == "CD",]
    }
    if(dataset_name == "Barlow") {
      counts_A <- data[groups == "control",]
      counts_B <- data[groups == "keto",]
    }
    if(dataset_name == "Song") {
      counts_A <- data[groups == "lung",]
      counts_B <- data[groups == "brain",]
    }
    if(dataset_name == "Monaco") {
      counts_A <- data[groups == "CD4_naive",]
      counts_B <- data[groups == "PBMC",]
    }
    if(dataset_name == "Hagai") {
      counts_A <- data[groups == "unstimulated",]
      counts_B <- data[groups == "pIC4",]
    }
    if(dataset_name == "Owens") {
      counts_A <- data[groups == "early_series",]
      counts_B <- data[groups == "late_series",]
    }
    if(dataset_name == "Klein") {
      counts_A <- data[groups == "unstimulated",]
      counts_B <- data[groups == "LIF-2hr",]
    }
    if(dataset_name == "Yu") {
      counts_A <- data[groups == "Brn",]
      counts_B <- data[groups == "Lvr",]
    }
  
    # This takes 2-3 min. to run on 15K features
    features <- as.data.frame(characterize_dataset(counts_A, counts_B))
  
    # Add a few more
    features$P <- ncol(counts_A)
  
    features <- features %>%
      select(c(-FW_RA_PFC1_D, FW_CLR_MED_D, FW_CLR_SD_D, FW_CLR_PNEG_D))
  
    real_features <- rbind(real_features,
                           features)
  }

  simulated_features <- pull_features()

  simulated_features <- simulated_features %>%
    filter(METHOD == "ALDEx2") %>%
    select(!c(METHOD, TPR, FPR))

  saveRDS(list(real_features = real_features,
               simulated_features = simulated_features),
          save_fn)
}

feature_intersection <- data.frame(feature = colnames(simulated_features),
                                   idx_sim = 1:ncol(simulated_features)) %>%
  inner_join(data.frame(feature = colnames(real_features),
                        idx_real = 1:ncol(real_features)), by = "feature")

mds_features <- rbind(simulated_features[,feature_intersection$idx_sim],
                      real_features[,feature_intersection$idx_real])
labels <- c(rep("simulated", nrow(simulated_features)),
            datasets)
sim_idx <- 1:nrow(simulated_features)
real_idx <- (nrow(simulated_features)+1):length(labels)

# Scale features
mds_features <- scale(mds_features)

# Testing
subsample_idx <- c(sample(sim_idx, size = 1000), real_idx)
mds_features <- mds_features[subsample_idx,]
labels <- labels[subsample_idx]

# 1) Diagnostic biplot

# features.pca <- prcomp(mds_features, scale. = FALSE)
# p0 <- ggbiplot(features.pca, obs.scale = 1, var.scale = 1,
#   groups = labels) +
#   scale_color_discrete(name = '')
# saveRDS(p0, "temp.rds")
# ggsave("temp.png",
#        p0,
#        dpi = 100,
#        units = "in",
#        height = 10,
#        width = 10)

# 2) Regular PCA plot

# Remove features related to totals
mds_features <- mds_features[,!(colnames(mds_features) %in% c(TOTALS_C_FC, TOTALS_C_D, TOTALS_C_MAX_D, TOTALS_C_MED_D, TOTALS_C_SD_D))]

feature_dist <- dist(mds_features)
coords <- cmdscale(feature_dist, k = 4)

plot_df <- data.frame(x1 = coords[,1], y1 = coords[,2],
                      x2 = coords[,3], y2 = coords[,4],
                      label = labels)
plot_df$label <- factor(plot_df$label,
                        levels = c("simulated", datasets))
levels(plot_df$label) <- c("simulated",
                           "Vieira-Silva et al.",
                           "Barlow et al.",
                           "Song et al.",
                           "Monaco et al.",
                           "Hagai et al.",
                           "Owens et al.",
                           "Klein et al.",
                           "Yu et al.")
p1 <- ggplot() +
  geom_point(data = plot_df %>% filter(label == "simulated"),
             mapping = aes(x = x1, y = y1),
             size = 3,
             color = "#cccccc") +
  geom_point(data = plot_df %>% filter(label != "simulated"),
             mapping = aes(x = x1, y = y1, fill = label),
             size = 3,
             shape = 21,
             stroke = 2) +
  scale_fill_ucscgb() +
  theme_bw() +
  labs(x = "PCA 1",
       y = "PCA 2",
       fill = "Data set")
ggsave(file.path("output", "images", paste0("features_PC1-2.png")),
       p1,
       dpi = 100,
       units = "in",
       height = 6,
       width = 7.5)

p2 <- ggplot() +
  geom_point(data = plot_df %>% filter(label == "simulated"),
             mapping = aes(x = x2, y = y2),
             size = 3,
             color = "#cccccc") +
  geom_point(data = plot_df %>% filter(label != "simulated"),
             mapping = aes(x = x2, y = y2, fill = label),
             size = 3,
             shape = 21,
             stroke = 2) +
  scale_fill_ucscgb() +
  theme_bw() +
  labs(x = "PCA 3",
       y = "PCA 4",
       fill = "Data set")
ggsave(file.path("output", "images", paste0("features_PC3-4.png")),
       p2,
       dpi = 100,
       units = "in",
       height = 6,
       width = 7.5)

difft <- Sys.time() - start
cat(paste0("Time elapsed: ", round(difft, 3), " ", attr(difft, "units"), "\n"))
