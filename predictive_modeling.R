source("path_fix.R")

library(codaDE)
library(tidyverse)
library(RSQLite)
library(mlegp)
library(tidyr)

output_dir <- file.path("output", "GP")
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}

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

for(use_method in c("NBGLM", "DESeq2", "ALDEx2", "MAST", "scran")) {
  for(use_result_type in c("fpr", "tpr")) {

    cat(paste0("Modeling ", use_result_type, " w/ ", use_method, "\n"))

    data <- results %>%
      filter(METHOD == use_method) %>%
      filter(RESULT_TYPE == use_result_type)

    # Subset for testing
    # data <- data[sample(1:nrow(data), size = 100),]

    # Scale the features
    uuids <- data$UUID
    features <- data %>%
      select(!c(UUID, METHOD, RESULT, RESULT_TYPE))
    features <- as.data.frame(apply(features, 2, function(x) {
      if(sd(x) > 0) {
        scale(x)
      } else {
        scale(x, scale = FALSE)
      }
    }))
    response <- data %>%
      select(RESULT)

    # Drop columns with no variation
    # In practice this only happens in testing/subsetting to small samples
    features <- features[,apply(features, 2, sd) > 0]

    # Define test/train set
    n <- nrow(data)
    train_idx <- sample(1:nrow(data), size = round(n*0.8))
    test_idx <- setdiff(1:nrow(data), train_idx)
    train_uuids <- uuids[train_idx]
    train_features <- features[train_idx,]
    train_response <- response[train_idx,]
    test_uuids <- uuids[test_idx]
    test_features <- features[test_idx,]
    test_response <- response[test_idx,]

    # Fit GP on training set
    start <- Sys.time()
    gp <- mlegp(train_features, train_response)
    cat(paste0("Elapsed fit time: ", Sys.time() - start, "\n"))

    # Predict on test set
    start <- Sys.time()
    output_pred <- predict(gp, newData = test_features, se.fit = FALSE)
    cat(paste0("Elapsed predict time: ", Sys.time() - start, "\n"))

    cat("Saving results...\n")
    saveRDS(list(train_features = cbind(UUID = train_uuids, train_features),
                 test_features = cbind(UUID = test_uuids, test_features),
                 train_response = train_response,
                 test_response = test_response,
                 prediction = output_pred[,1]),
            file.path(output_dir,
                      paste0("results_", use_method, "_", use_result_type, ".rds")))
    
    plot_df <- data.frame(true = test_response,
                          predicted = output_pred[,1],
                          p = factor(test_features$P))
    levels(plot_df$p) <- c("low", "med", "high")
    pl <- ggplot(plot_df, aes(x = true, y = predicted, fill = p)) +
      geom_point(shape = 21, size = 3) +
      labs(x = paste0("observed ", toupper(use_result_type)),
           y = paste0("predicted ", toupper(use_result_type)),
           fill = "feature no.")
    ggsave(file.path("output", "images", paste0("GP_predictions_", use_method, "_", use_result_type, ".png")),
           plot = pl,
           dpi = 100,
           units = "in",
           height = 6,
           width = 7)

    dbDisconnect(conn)
  }
}
