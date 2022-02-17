source("path_fix.R")

library(tidyverse)
library(codaDE)
library(ggsci)
library(RColorBrewer)
# library(ggbiplot)
# library("devtools")
# install_github("vqv/ggbiplot")

all_features <- readRDS(file.path("output", "feature_dump_cpm.rds"))

pull_top_features <- function(method) {
  top_features_fpr <- readRDS(file.path("output",
                                        "predictive_fits",
                                        "oracle",
                                        "regression_cpm",
                                        paste0("oracle_FPR_",method,"_variable-importance.rds")))
  top_features_fpr$feature <- rownames(top_features_fpr)
  top_features_fpr <- top_features_fpr %>%
    arrange(desc(Overall)) %>%
    slice(1:5) %>%
    pull(feature)
  top_features_tpr <- readRDS(file.path("output",
                                        "predictive_fits",
                                        "oracle",
                                        "regression_cpm",
                                        paste0("oracle_TPR_",method,"_variable-importance.rds")))
  top_features_tpr$feature <- rownames(top_features_tpr)
  top_features_tpr <- top_features_tpr %>%
    arrange(desc(Overall)) %>%
    slice(1:5) %>%
    pull(feature)
  sort(unique(c(top_features_tpr, top_features_fpr)))
}

top_features <- unique(c(pull_top_features("ALDEx2"), pull_top_features("DESeq2"), pull_top_features("scran")))
mds_features <- all_features %>%
  select(top_features) %>%
  mutate(type = "simulation")
# Subsample for testing
mds_features <- mds_features[sample(1:nrow(mds_features), size = 2000),]

datasets <- c("Hagai", "Monaco", "Song", "Hashimshony", "Barlow", "Gruen",
              "Muraro", "Owens", "VieiraSilva", "Kimmerling", "Yu", "Klein")
thresholds <- rep(1, length(datasets))

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
coords <- cmdscale(feature_dist, k = 4) # Takes about 10-20 sec.

# 1) Plot simulations in space of top 5-ish features
plot_df <- data.frame(x = coords[,1], y = coords[,2], z1 = coords[,3], z2 = coords[,4], type = mds_features$type)
plot_df$type <- factor(plot_df$type, levels = c(datasets, "simulation"))

palette <- generate_highcontrast_palette(length(datasets))

p1 <- ggplot() +
  geom_boxplot(data = mds_features %>% pivot_longer(!type, names_to = "feature", values_to = "value") %>% filter(type == "simulation"),
               mapping = aes(x = 1, y = value), outlier.shape = NA) +
  geom_point(data = mds_features %>% pivot_longer(!type, names_to = "feature", values_to = "value") %>% filter(type != "simulation"),
             mapping = aes(x = 1, y = value, fill = type), size = 2, shape = 21) +
  facet_wrap(. ~ feature, scales = "free_y") +
  scale_fill_manual(values = palette) +
  theme_bw()

show(p1)

ggsave(file.path("output", "images", paste0("MDS_axes_1-2.png")),
       p1,
       dpi = 100,
       units = "in",
       height = 5,
       width = 6.5)
