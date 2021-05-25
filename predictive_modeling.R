source("path_fix.R")

library(codaDE)
library(tidyverse)
library(RSQLite)
library(mlegp)
library(tidyr)

# Create DB (if doesn't exist?)
conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

results <- dbGetQuery(conn, paste0("SELECT datasets.UUID, P, PARTIAL, ",
                                   "FOLD_CHANGE, MEAN_CORR, MEDIAN_CORR, ",
                                   "BASE_SPARSITY, DELTA_SPARSITY, ",
                                   "PERCENT_STABLE, DIR_CONSENSUS, ",
                                   "MAX_DELTA_REL, MEDIAN_DELTA_REL, ",
                                   "BASE_ENTROPY, DELTA_ENTROPY, METHOD, ",
                                   "RESULT, RESULT_TYPE FROM ",
                                   "datasets LEFT JOIN characteristics ",
                                   "ON datasets.UUID=characteristics.UUID ",
                                   "LEFT JOIN results ",
                                   "ON (characteristics.UUID=results.UUID ",
                                   "AND characteristics.partial=results.PARTIAL_INFO);"))

use_method <- "NBGLM"
use_result_type <- "tpr"

data <- results %>%
  filter(METHOD == use_method) %>%
  filter(RESULT_TYPE == use_result_type)


# Scale the features
features <- data %>%
  select(!c(UUID, METHOD, RESULT, RESULT_TYPE))
features <- apply(features, 2, function(x) {
  if(sd(x) > 0) {
    scale(x)
  } else {
    scale(x, scale = FALSE)
  }
})
response <- data %>%
  select(RESULT)

# Drop columns with no variation
# In practice this only happens in testing/subsetting to small samples
features <- features[,apply(features, 2, sd) > 0]

# Define test/train set
n <- nrow(data)
train_idx <- sample(1:nrow(data), size = round(n*0.8))
test_idx <- setdiff(1:nrow(data), train_idx)
train_features <- features[train_idx,]
train_response <- response[train_idx,]
test_features <- features[test_idx,]
test_response <- response[test_idx,]

# Fit GP on training set
start <- Sys.time()
gp <- mlegp(train_features, train_response)
cat(paste0("Elapsed time: ", Sys.time() - start, "\n"))

# Predict on test set
start <- Sys.time()
output_pred <- predict(gp, newData = test_features, se.fit = FALSE)
cat(paste0("Elapsed time: ", Sys.time() - start, "\n"))

plot_df <- data.frame(true = test_response,
                      predicted = output_pred[,1])
ggplot(plot_df, aes(x = true, y = predicted)) +
  geom_point(shape = 21, size = 3) +
  labs(x = "observed FPR",
       y = "predicted FPR")

# Evaluate on test (visualize)

dbDisconnect(conn)
