source("path_fix.R")

library(codaDE)
library(RSQLite)
library(mlegp)
library(tidyr)

# Create DB (if doesn't exist?)
conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

# Pull data characteristics and results (join)
data <- dbGetQuery(conn, paste0("SELECT P, PARTIAL_INFO, ",
                                "FOLD_CHANGE, MEAN_CORR, ",
                                "MEDIAN_CORR, RESULT ",
                                "FROM results LEFT JOIN datasets ",
                                "ON results.UUID=datasets.UUID ",
                                "WHERE RESULT_TYPE='fpr' AND METHOD='NBGLM' ",
                                "AND PARTIAL_INFO=0;"))

data <- bind_rows(data,
                  dbGetQuery(conn, paste0("SELECT P, PARTIAL_INFO, ",
                                          "FOLD_CHANGE_PARTIAL AS FOLD_CHANGE, ",
                                          "MEAN_CORR_PARTIAL AS MEAN_CORR, ",
                                          "MEDIAN_CORR_PARTIAL AS MEDIAN_CORR, ",
                                          "RESULT FROM ",
                                          "results LEFT JOIN datasets ",
                                          "ON results.UUID=datasets.UUID ",
                                          "WHERE RESULT_TYPE='fpr' ",
                                          "AND METHOD='NBGLM' ",
                                          "AND PARTIAL_INFO=1;")))
partial_info_labels <- data$PARTIAL_INFO
p_labels <- data$P
data <- data %>%
  select(!PARTIAL_INFO)
head(data)

# Scale the features
features <- data %>%
  select(!RESULT)
features <- apply(features, 2, scale)
response <- data %>%
  select(RESULT)

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

plot_df <- data.frame(true = test_set$RESULT,
                      predicted = output_pred[,1],
                      # label = factor(partial_info_labels[test_idx]))
                      label = factor(p_labels[test_idx]))
ggplot(plot_df, aes(x = true, y = predicted, fill = label)) +
  geom_point(shape = 21, size = 3) +
  labs(x = "observed FPR",
       y = "predicted FPR")

# Evaluate on test (visualize)

dbDisconnect(conn)
