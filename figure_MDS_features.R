source("path_fix.R")

library(tidyverse)
library(codaDE)
library(ggsci)
library(ggbiplot)
library(RColorBrewer)
# library("devtools")
# install_github("vqv/ggbiplot")

all_features <- readRDS(file.path("output", "feature_dump.rds"))

pull_top_features <- function(method) {
  top_features <- readRDS(file.path("output",
                                    "predictive_fits",
                                    "oracle",
                                    "regression_baseline",
                                    paste0("oracle_FPR_",method,"_variable-importance.rds")))
  top_features$feature <- rownames(top_features)
  top_features <- top_features %>%
    arrange(desc(Overall)) %>%
    slice(1:5) %>%
    pull(feature)
  top_features
}

save_fn <- file.path("output", "MDS_features.rds")
if(!file.exists(save_fn)) {
  top_features <- unique(c(pull_top_features("ALDEx2"), pull_top_features("DESeq2"), pull_top_features("scran")))
  mds_features <- all_features %>%
    select(top_features) %>%
    mutate(type = "simulation")
  # mds_features <- mds_features[sample(1:nrow(mds_features), size = 1000),] # subsample for testing
  
  datasets <- c("VieiraSilva", "Barlow", "Song", "Monaco", "Hagai", "Owens", "Klein", "Yu")
  thresholds <- c(1, 1, 1, 1, 1, 1, 1, 1)
  
  for(i in 1:length(datasets)) {
    this_dataset <- datasets[i]
    this_threshold <- thresholds[i]
    temp <- readRDS(file.path("output",
                              paste0("filtered_features_",this_dataset,"_threshold",this_threshold,".rds"))) %>%
      select(top_features) %>%
      mutate(type = this_dataset)
    mds_features <- rbind(mds_features, temp)
  }
  
  feature_dist <- dist(mds_features[,1:(ncol(mds_features)-1)])
  coords <- cmdscale(feature_dist, k = 4)
  
  # 1) Plot simulations in space of top 5-ish features
  plot_df <- data.frame(x = coords[,1], y = coords[,2], z1 = coords[,3], z2 = coords[,4], type = mds_features$type)
  plot_df$type <- factor(plot_df$type, levels = c("Barlow",
                                                  "VieiraSilva",
                                                  "Song",
                                                  "Hagai",
                                                  "Monaco",
                                                  "Owens",
                                                  "Klein",
                                                  "Yu",
                                                  "simulation"))
  levels(plot_df$type) <- c("Barlow et al.",
                            "Vieira-Silva et al.",
                            "Song et al.",
                            "Hagai et al.",
                            "Monaco et al.",
                            "Owens et al.",
                            "Klein et al.",
                            "Yu et al.",
                            "simulation")
  
  saveRDS(plot_df, save_fn)
} else {
  plot_df <- readRDS(save_fn)
}

p1 <- ggplot() +
  geom_point(data = plot_df %>% filter(type == "simulation"),
             mapping = aes(x = x, y = y),
             color = "gray") +
  geom_point(data = plot_df %>% filter(type != "simulation"),
             mapping = aes(x = x, y = y, fill = type),
             size = 3,
             shape = 21) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(x = "PC 1",
       y = "PC 2",
       fill = "Data set")
# p2 <- ggplot() +
#   geom_point(data = plot_df %>% filter(type == "simulation"),
#              mapping = aes(x = z1, y = z2),
#              color = "gray") +
#   geom_point(data = plot_df %>% filter(type != "simulation"),
#              mapping = aes(x = z1, y = z2, fill = type),
#              size = 3,
#              shape = 21) +
#   scale_fill_brewer(palette = "Set1") +
#   theme_bw() +
#   labs(x = "PC 1",
#        y = "PC 2",
#        fill = "Data set")

ggsave(file.path("output", "images", paste0("features_PC1-2.png")),
       p1,
       dpi = 100,
       units = "in",
       height = 5,
       width = 6.5)
