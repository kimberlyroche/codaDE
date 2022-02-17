library(codaDE)
library(tidyr)
library(ggplot2)
library(scales)
library(dplyr)
library(glmnet)
library(caret)
library(doParallel)
registerDoParallel(4)

plot_feature_corr <- function(local_features) {
  feature_corr <- cor(local_features, method = "spearman")
  feature_corr <- cbind(rownames(feature_corr), feature_corr)
  colnames(feature_corr) <- c("feature1", rownames(feature_corr))
  rownames(feature_corr) <- NULL
  feature_corr <- as.data.frame(feature_corr)
  feature_corr <- pivot_longer(feature_corr,
                               !feature1,
                               names_to = "feature2",
                               values_to = "value")
  feature_corr$value <- as.numeric(feature_corr$value)
  ggplot(feature_corr, aes(x = feature1, y = feature2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = muted("blue"),
                         mid = "white",
                         high = muted("red"),
                         midpoint = 0) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(fill = "Spearman\nrho")
}

abs_feature_list <- c("LOG_MEAN", "PERTURBATION",
                      "REP_NOISE", "FC_ABSOLUTE", "MED_ABS_TOTAL",
                      "MED_REL_TOTAL", "PERCENT_DIFF_SIM", "PERCENT_DIFF_REALIZ")

features <- pull_features(DE_methods = c("ALDEx2", "DESeq2", "NBGLM", "scran"),
                          use_baseline = "self",
                          exclude_partials = TRUE,
                          exclude_independent = FALSE,
                          feature_list = NULL,
                          abs_feature_list = abs_feature_list)

features <- features %>%
  select(!TPR)

# ------------------------------------------------------------------------------
#   How do FC and PERCENT DE associate with FPR?
# ------------------------------------------------------------------------------

temp <- features %>%
  filter(METHOD == "scran") %>%
  select(PERCENT_DIFF_REALIZ, REP_NOISE, FC_ABSOLUTE, FPR, P)
temp$FC_ABSOLUTE[temp$FC_ABSOLUTE < 1] <- 1/temp$FC_ABSOLUTE[temp$FC_ABSOLUTE < 1]
temp <- temp %>%
  filter(FC_ABSOLUTE < 5)

# Visualize percent DE
ggplot(data.frame(x = temp$FC_ABSOLUTE, y = temp$FPR, p = temp$P, label = temp$PERCENT_DIFF_REALIZ),
       aes(x = x, y = y, fill = label)) +
  geom_point(size = 2, shape = 21) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = mean(temp$PERCENT_DIFF_REALIZ)) +
  facet_wrap(. ~ p) +
  theme_bw() +
  labs(x = "fold change (in direction of increase)",
       y = "FPR",
       fill = "% diff. features")

# ------------------------------------------------------------------------------
#   ID data sets where we have (1) high FPR (2) <25% DE (3) low FC
# ------------------------------------------------------------------------------

temp <- features %>%
  filter(METHOD == "NBGLM") %>%
  filter(FPR > 0.5) %>%
  filter(PERCENT_DIFF_REALIZ < 0.25) %>%
  filter(REP_NOISE < 0.5)
temp$FC_ABSOLUTE[temp$FC_ABSOLUTE < 1] <- 1/temp$FC_ABSOLUTE[temp$FC_ABSOLUTE < 1]
temp <- temp %>%
  filter(FC_ABSOLUTE < 1.5)

idx1 <- sample(1:nrow(temp), size = 1)
uuid <- temp$UUID[idx1]
temp$P[idx1]
temp$FC_ABSOLUTE[idx1]
uuid
data <- readRDS(file.path("output", "datasets", paste0(uuid, ".rds")))
palette <- generate_highcontrast_palette(1000)
plot_stacked_bars(data$simulation$abundances, palette)
plot_stacked_bars(data$simulation$observed_counts1, palette)

# Pull the calls for this data set
conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
baseline_calls <- dbGetQuery(conn, paste0("SELECT BASELINE_CALLS FROM datasets WHERE UUID='", uuid, "'"))$BASELINE_CALLS
calls <- dbGetQuery(conn, paste0("SELECT CALLS FROM results WHERE UUID='", uuid, "' AND METHOD='NBGLM'"))$CALLS
dbDisconnect(conn)

baseline_calls <- as.numeric(strsplit(baseline_calls, ";")[[1]])
calls <- as.numeric(strsplit(calls, ";")[[1]])

alpha <- 0.05 / length(baseline_calls)
tp <- which(baseline_calls < alpha & calls < alpha)
fp <- which(baseline_calls > alpha & calls < alpha)

idx <- sample(fp, size = 1)
plot_df <- data.frame(x = rep(1:nrow(data$simulation$abundances), 2),
                  y = c(data$simulation$abundances[,idx], data$simulation$observed_counts1[,idx]),
                  label1 = factor(rep(data$simulation$groups, 2)),
                  label2 = c(rep("absolute", nrow(data$simulation$abundances)),
                             rep("observed", nrow(data$simulation$abundances))))
levels(plot_df$label1) <- c("A", "B")
ggplot(plot_df, aes(x = x, y = y, fill = label1)) +
  geom_point(size = 3, shape = 21) +
  facet_wrap(. ~ label2, scales = "free_y") +
  theme_bw() +
  labs(x = "sample index",
       y = "abundance",
       fill = "Group")

# ------------------------------------------------------------------------------

abs_features <- features %>%
  select(c(abs_feature_list, "P", "METHOD", "FPR"))

rel_features <- features %>%
  select(!abs_feature_list)

# p <- plot_feature_corr(abs_features[,1:9])
# p <- plot_feature_corr(rel_features[,1:57])
p <- plot_feature_corr(features[,1:65])
show(p)
ggsave(file.path("output", "images", "correlated_features.png"),
       dpi = 100,
       units = "in",
       height = 9,
       width = 10)

# ------------------------------------------------------------------------------
#   LM and LR models on relative abundance data
#
#   This is a great reference:
#   http://www.sthda.com/english/articles/36-classification-methods-essentials/149-penalized-logistic-regression-essentials-in-r-ridge-lasso-and-elastic-net/
# ------------------------------------------------------------------------------

feature_df <- abs_features
feature_df <- rel_features
feature_df <- features

# Remove features apparently correlated with absolute information
# feature_df <- abs_features
# feature_df <- rel_features
feature_df <- features
feature_df <- feature_df %>%
  select(!c(FW_RA_PFC05_D, FW_CLR_PFC05_D))
feature_df$FC_ABSOLUTE[feature_df$FC_ABSOLUTE < 1] <- 1/feature_df$FC_ABSOLUTE[feature_df$FC_ABSOLUTE < 1]
# Make binary
feature_df$FC_ABSOLUTE <- sapply(feature_df$FC_ABSOLUTE, function(x) {
  ifelse(x > 2, 1, 0)
})

for(method_of_choice in c("ALDEx2", "DESeq2", "scran")) {
# for(method_of_choice in c("NBGLM")) {
  # Per-model for now
  if(method_of_choice == "all") {
    method <- feature_df$METHOD
    lr_features <- feature_df %>%
      select(!c(METHOD))
  } else {
    lr_features <- feature_df %>%
      filter(METHOD == method_of_choice) %>%
      select(!c(METHOD))
  }
  
  lr_feat_scaled <- lr_features
  lr_feat_scaled[,1:(ncol(lr_features)-1)][,1:(ncol(lr_features)-1)] <- apply(lr_feat_scaled[,1:(ncol(lr_features)-1)], 2, scale)
  
  if(method_of_choice == "all") {
    lr_feat_scaled$METHOD <- method
  }
  
  x <- model.matrix(FPR ~ ., lr_feat_scaled)[,-1]
  y1 <- lr_feat_scaled$FPR # regression outcomes
  y2 <- ifelse(lr_feat_scaled$FPR > 0.05, 1, 0) # classification outcomes
  
  # LM fit
  cvfit <- cv.glmnet(x, y1,
                     family = "gaussian",
                     type.measure = "mse",
                     alpha = 1,
                     parallel = TRUE)
  lambda_idx <- which(cvfit$lambda == cvfit$lambda.1se)
  cat(paste0(method_of_choice, "\n"))
  cat(paste0("\tLM # features (best lambda): ", cvfit$nzero[lambda_idx], "\n"))
  
  # Gauge accuracy as R^2 on predictions
  train_idx <- createDataPartition(y1,
                                   p = 0.8,
                                   list = FALSE,
                                   times = 100)
  r2_vec <- c()
  for(j in 1:ncol(train_idx)) {
    fit <- glmnet(x[train_idx[,j],],
                  y1[train_idx[,j]],
                  family = "gaussian",
                  alpha = 1,
                  lambda = cvfit$lambda.1se)
    predictions <- predict(fit, x[-train_idx[,j],])
  
    r2 <- cor(y1[-train_idx[,j]], as.vector(predictions))**2
    r2_vec <- c(r2_vec, r2)
  }
  cat(paste0("\tLM R^2 (best lambda): ", round(mean(r2_vec), 2), "\n"))
  
  # Plot combined features retained
  # betas <- as.matrix(coef(fit))
  # retained_idx <- which(betas != 0)
  # retained_idx <- retained_idx[2:length(retained_idx)]
  # retained_idx <- retained_idx - 1
  # colnames(x)[retained_idx]

  # p <- plot_feature_corr(x[,retained_idx])
  # show(p)
  # ggsave(file.path("output", "images", "correlated_features_DESeq2.png"),
  #        dpi = 100,
  #        units = "in",
  #        height = 6,
  #        width = 7)
  
  # LR fit
  # cvfit <- cv.glmnet(x, y2,
  #                    family = "binomial",
  #                    type.measure = "class",
  #                    alpha = 1,
  #                    parallel = TRUE)
  # lambda_idx <- which(cvfit$lambda == cvfit$lambda.1se)
  # cat(paste0("\tLR # features (best lambda): ", cvfit$nzero[lambda_idx], "\n"))
  # cat(paste0("\tLR classification accuracy (best lambda): ", round(1 - cvfit$cvm[lambda_idx], 2), "\n"))
}

# ------------------------------------------------------------------------------
#   Look at correlated features
# ------------------------------------------------------------------------------

plot_df <- data.frame(x = features$FC_ABSOLUTE, y = features$FW_RA_PFC05_D)
ggplot(plot_df, aes(x = x, y = y)) +
  geom_point(size = 2, shape = 21, fill = "#bbbbbb") +
  theme_bw() +
  xlim(0, 100) +
  labs(x = "absolute fold change",
       y = "percent features w/ < 0.5 FC in rel. abundance",
       title = paste0("Spearman rho = ", round(cor(plot_df$x, plot_df$y, method = "spearman"), 3)))


plot_df <- data.frame(x = features$PERCENT_DIFF_REALIZ, y = features$FW_CLR_PFC05_D)
ggplot(plot_df, aes(x = x, y = y)) +
  geom_point(size = 2, shape = 21, fill = "#bbbbbb") +
  theme_bw() +
  labs(x = "percent differentially abundant features (by NB GLM)",
       y = "percent features with < 0.5 FC in CLR",
       title = paste0("Spearman rho = ", round(cor(plot_df$x, plot_df$y, method = "spearman"), 3)))








