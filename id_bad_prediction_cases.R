library(codaDE)
library(tidyverse)

plot_labeling <- function(plot_df, label_idx) {
  plot_df$selected <- FALSE
  plot_df$selected[label_idx] <- TRUE
  ggplot(plot_df, aes(x = x, y = y, fill = selected)) +
    geom_point(size = 3, shape = 21)
}

data <- readRDS("output/GP/results_DESeq2_fpr.rds")
palette <- generate_highcontrast_palette(100)

# Pull low P cases
low_p_feat_idx <- which(data$test_features$P < -0.8)
plot_df <- data.frame(x = data$test_response[low_p_feat_idx],
                      y = data$prediction[low_p_feat_idx])

# The poorly predicted low-P cases are both PARTIAL = 0 and PARTIAL = 1
# plot_labeling(plot_df, intersect(low_p_feat_idx,
#                                  which(data$test_features$PARTIAL > 0)))

# Pull "poorly predicted" (flat-line) cases
delta <- abs(0.271 - data$prediction[low_p_feat_idx])
poorly_predicted <- which(delta < 0.0025)
plot_labeling(plot_df, poorly_predicted)

# Pull especially bad predictions (subset of these)
bad_ref <- which(delta < 0.0025 &
                   (data$test_response[low_p_feat_idx] < 0.2 |
                      data$test_response[low_p_feat_idx] > 0.4))

plot_labeling(plot_df, bad_ref)

# Pull especially good predictions
delta <- abs(data$test_response[low_p_feat_idx] - data$prediction[low_p_feat_idx])
good_ref <- which(delta < 0.1) # not "well" but ok
good_ref <- setdiff(good_ref, poorly_predicted)

plot_labeling(plot_df, good_ref)

# Get UUIDs
uuid_bad <- data$test_features$UUID[bad_ref]
uuid_good <- data$test_features$UUID[well_predicted]

# Visualize a "bad" case
idx <- sample(1:length(uuid_bad), size = 1)
bad_data <- readRDS(file.path("output", "datasets", paste0(uuid_bad[idx], ".rds")))
plot_stacked_bars(bad_data$simulation$abundances, palette = palette)
plot_stacked_bars(bad_data$simulation$observed_counts1, palette = palette)
data$test_features %>% filter(UUID == uuid_bad[idx])

# Visualize a "good" case
idx <- sample(1:length(uuid_good), size = 1)
bad_data <- readRDS(file.path("output", "datasets", paste0(uuid_good[idx], ".rds")))
# plot_stacked_bars(bad_data$simulation$abundances, palette = palette)
plot_stacked_bars(bad_data$simulation$observed_counts1, palette = palette)

# Variation in features in the bad cases vs. good
features_bad <- data$test_features %>% filter(UUID %in% uuid_bad)
features_good <- data$test_features %>% filter(UUID %in% uuid_good)

# feature_name <- sample(colnames(data$test_features)[4:ncol(data$test_features)],
#                        size = 1)
feature_name <- "DIR_CONSENSUS"
plot_df <- data.frame(feature = c(features_bad[[feature_name]],
                                  features_good[[feature_name]]),
                      type = c(rep("bad", nrow(features_bad)),
                               rep("good", nrow(features_good))))
ggplot(plot_df, aes(x = feature, fill = type)) +
  geom_histogram(alpha = 1, color = "black") +
  facet_wrap(. ~ type) +
  labs(title = feature_name)


