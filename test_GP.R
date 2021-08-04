source("path_fix.R")

library(codaDE)
library(RSQLite)
library(tidyverse)
library(mlegp)

# GP models take a while to run. Save the results for later tinkering.
save_path <- list("output", "predictive_fits", "all", "GP_results")
for(i in 1:length(save_path)) {
  fp <- do.call(file.path, save_path[1:i])
  if(!dir.exists(fp)) {
    dir.create(fp)
  }
}
save_dir <- fp

use_baseline <- "self"
use_result_type <- "fpr"
testing <- TRUE

save_fn <- "temp.rds"
if(file.exists(save_fn)) {
  data <- readRDS(save_fn)
} else {
  conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
  data <- dbGetQuery(conn, paste0("SELECT datasets.UUID, P, CORRP, ",
                                  "TOTALS_C_FC, ",
                                  "TOTALS_C_D, ",
                                  "TOTALS_C_MAX_D, ",
                                  "TOTALS_C_MED_D, ",
                                  "TOTALS_C_SD_D, ",
                                  "CORR_RA_MED, ",
                                  "CORR_RA_SD, ",
                                  "CORR_RA_SKEW, ",
                                  "CORR_LOG_MED, ",
                                  "CORR_LOG_SD, ",
                                  "CORR_LOG_SKEW, ",
                                  "CORR_CLR_MED, ",
                                  "CORR_CLR_SD, ",
                                  "CORR_CLR_SKEW, ",
                                  "COMP_C_P0_A, ",
                                  "COMP_C_P0_B, ",
                                  "COMP_C_P1_A, ",
                                  "COMP_C_P1_B, ",
                                  "COMP_C_P5_A, ",
                                  "COMP_C_P5_B, ",
                                  "COMP_RA_P01_A, ",
                                  "COMP_RA_P01_B, ",
                                  "COMP_RA_P1_A, ",
                                  "COMP_RA_P1_B, ",
                                  "COMP_RA_P5_A, ",
                                  "COMP_RA_P5_B, ",
                                  "COMP_RA_MAX_A, ",
                                  "COMP_RA_MED_A, ",
                                  "COMP_RA_SD_A, ",
                                  "COMP_RA_SKEW_A, ",
                                  "COMP_RA_MAX_B, ",
                                  "COMP_RA_MED_B, ",
                                  "COMP_RA_SD_B, ",
                                  "COMP_RA_SKEW_B, ",
                                  "COMP_C_ENT_A, ",
                                  "COMP_C_ENT_B, ",
                                  "FW_RA_MAX_D, ",
                                  "FW_RA_MED_D, ",
                                  "FW_RA_SD_D, ",
                                  "FW_RA_PPOS_D, ",
                                  "FW_RA_PNEG_D, ",
                                  "FW_RA_PFC05_D, ",
                                  "FW_RA_PFC1_D, ",
                                  "FW_RA_PFC2_D, ",
                                  "FW_LOG_MAX_D, ",
                                  "FW_LOG_MED_D, ",
                                  "FW_LOG_SD_D, ",
                                  "FW_LOG_PPOS_D, ",
                                  "FW_LOG_PNEG_D, ",
                                  "FW_LOG_PFC05_D, ",
                                  "FW_LOG_PFC1_D, ",
                                  "FW_LOG_PFC2_D, ",
                                  "FW_CLR_MAX_D, ",
                                  "FW_CLR_MED_D, ",
                                  "FW_CLR_SD_D, ",
                                  "FW_CLR_PPOS_D, ",
                                  "FW_CLR_PNEG_D, ",
                                  "FW_CLR_PFC05_D, ",
                                  "FW_CLR_PFC1_D, ",
                                  "FW_CLR_PFC2_D, ",
                                  "METHOD, TPR, FPR FROM ",
                                  "datasets LEFT JOIN characteristics ",
                                  "ON datasets.UUID=characteristics.UUID ",
                                  "LEFT JOIN results ",
                                  "ON (characteristics.UUID=results.UUID ",
                                  "AND characteristics.partial=results.PARTIAL_INFO) ",
                                  "WHERE results.PARTIAL_INFO=0 AND ",
                                  "BASELINE_TYPE='", use_baseline, "'",
                                  ifelse(testing, " LIMIT 100", ""), ";"))
  dbDisconnect(conn)
  saveRDS(data, save_fn)
}

# Drop TPR or FPR missing results. These are typically funky simulations with
# effectively no differential expression.
data <- data %>%
  filter(!is.na(TPR) & !is.na(FPR))

# These predictors appear to be strongly correlated.
# Let's drop them for now.
data <- data %>%
  select(-c(FW_RA_PFC1_D, FW_CLR_MED_D, FW_CLR_SD_D, FW_CLR_PNEG_D))

# Scale the features
features <- data %>%
  select(!c(UUID, TPR, FPR))
features$METHOD <- factor(features$METHOD)

factors <- which(colnames(features) %in% c("METHOD"))
non_factors <- setdiff(1:ncol(features), factors)
features_nonfactors <- features[,non_factors]
features_nonfactors <- as.data.frame(apply(features_nonfactors, 2, function(x) {
  if(sd(x) > 0) {
    scale(x)
  } else {
    scale(x, scale = FALSE)
  }
}))
# Add fake variance to invariant features; just for testing, where subsets of
# samples are so small otherwise informative features (e.g. P) don't vary.
if(testing) {
  invariant_idx <- which(apply(features_nonfactors, 2, sd) == 0)
  for(i in invariant_idx) {
    features_nonfactors[,i] <- features_nonfactors[,i] +
      rbinom(nrow(features_nonfactors), size = 1, prob = 0.2)
  }
}

features <- cbind(METHOD = features[,factors], features_nonfactors)

write.table(features, "features.txt")
write.table(data %>% select(TPR), "response_TPR.txt")
write.table(data %>% select(FPR), "response_FPR.txt")

quit()

if(use_result_type == "tpr") {
  response <- data %>%
    select(TPR)
} else {
  response <- data %>%
    select(FPR)
}

n <- nrow(data)

uuids <- data$UUID

# Define test/train set for all remaining methods
train_idx <- sample(1:nrow(data), size = round(n*0.8))
test_idx <- setdiff(1:nrow(data), train_idx)
train_uuids <- uuids[train_idx]
train_features <- features[train_idx,]
train_response <- response[train_idx,]
if(use_result_type == "fpr") {
  train_response <- 1 - train_response
}
test_uuids <- uuids[test_idx]
test_features <- features[test_idx,]
test_response <- response[test_idx,]
if(use_result_type == "fpr") {
  test_response <- 1 - test_response
}

save_fn <- file.path(save_dir,
                     paste0("all_",
                            use_result_type,
                            "_",
                            use_baseline,
                            ".rds"))
if(!file.exists(save_fn)) {
  start <- Sys.time()
  train_features$METHOD <- as.numeric(train_features$METHOD)
  res <- mlegp(train_features, train_response)
  cat(paste0("Elapsed fit time: ", Sys.time() - start, "\n"))
  saveRDS(list(result = res,
               train_features = train_features,
               train_response = train_response,
               test_features = test_features,
               test_response = test_response), save_fn)
} else {
  res_obj <- readRDS(save_fn)
  res <- res_obj$result
  train_features <- res_obj$train_features
  train_response <- res_obj$train_response
  test_features <- res_obj$test_features
  test_response <- res_obj$test_response
}

start <- Sys.time()
test_features$METHOD <- as.numeric(test_features$METHOD)
output_pred <- predict(res, newData = test_features, se.fit = FALSE)
cat(paste0("Elapsed predict time: ", Sys.time() - start, "\n"))

# Render plot
plot_df <- data.frame(true = test_response,
                      predicted = output_pred[,1],
                      p = test_features$P)

plot_labels <- list(fpr = "specificity (1 - FPR)", tpr = "sensitivity (TPR)")

pl <- ggplot(plot_df, aes(x = true, y = predicted, fill = factor(p))) +
  geom_point(shape = 21, size = 3) +
  scale_fill_brewer(palette = "RdYlGn") +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  labs(x = paste0("observed ", plot_labels[[use_result_type]]),
       y = paste0("predicted ", plot_labels[[use_result_type]]),
       fill = "Feature number")
ggsave(file.path(save_dir,
                 paste0("predictions_all_",
                        use_result_type,
                        "_",
                        use_baseline,
                        ".png")),
       plot = pl,
       dpi = 100,
       units = "in",
       height = 4,
       width = 5.5)

# Calculate R^2
cat(paste0("R^2 (",toupper(use_result_type),"): ",
           round(cor(plot_df$true, plot_df$predicted)^2, 3), "\n"))
