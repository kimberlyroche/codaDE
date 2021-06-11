source("path_fix.R")

library(codaDE)
library(tidyverse)
library(RSQLite)
library(mlegp)
library(randomForest)
library(tidyr)
library(optparse)

# For plotting betas
library(jtools)
library(ggstance)
library(broom.mixed)

option_list = list(
  make_option(c("--model"), type = "character", default = "GP",
              help = "predictive model to use: GP, linear, RF", metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);
model <- opt$model
if(!(model %in% c("GP", "linear", "RF"))) {
  stop("Invalid model specified!\n")
}

# GP models take a while to run. Save the results for later tinkering.
output_dir <- "output"
if(model == "GP") {
  output_dir <- file.path("output", model)
  if(!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
}

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

mean_RMSE <- 0

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

    if(model == "linear") {
      lm_train_data <- cbind(train_features, train_response)
      res <- lm(train_response ~ P + PARTIAL + FOLD_CHANGE + MEAN_CORR + MEDIAN_CORR + BASE_SPARSITY + DELTA_SPARSITY + PERCENT_STABLE + DIR_CONSENSUS + MAX_DELTA_REL + MEDIAN_DELTA_REL + BASE_ENTROPY + DELTA_ENTROPY, data = lm_train_data)
      plot_summs(res)
      lm_test_data <- cbind(test_features, test_response)
      res_pred <- predict(res, newdata = lm_test_data)
      plot_df <- data.frame(true = test_response,
                            predicted = res_pred,
                            p = factor(test_features$P))
      levels(plot_df$p) <- c("low", "med", "high")
      pl <- ggplot(plot_df, aes(x = true, y = predicted, fill = p)) +
        geom_point(shape = 21, size = 3) +
        labs(x = paste0("observed ", toupper(use_result_type)),
             y = paste0("predicted ", toupper(use_result_type)),
             fill = "feature no.")
    }

    if(model == "RF") {
      rf_train_data <- cbind(train_features, train_response)
      res <- randomForest(train_response ~ P + PARTIAL + FOLD_CHANGE + MEAN_CORR + MEDIAN_CORR + BASE_SPARSITY + DELTA_SPARSITY + PERCENT_STABLE + DIR_CONSENSUS + MAX_DELTA_REL + MEDIAN_DELTA_REL + BASE_ENTROPY + DELTA_ENTROPY, data = rf_train_data)
      rf_test_data <- cbind(test_features, test_response)
      res_pred <- predict(res, newdata = rf_test_data)
      plot_df <- data.frame(true = test_response,
                            predicted = res_pred,
                            p = factor(test_features$P))
      levels(plot_df$p) <- c("low", "med", "high")
      pl <- ggplot(plot_df, aes(x = true, y = predicted, fill = p)) +
        geom_point(shape = 21, size = 3) +
        labs(x = paste0("observed ", toupper(use_result_type)),
             y = paste0("predicted ", toupper(use_result_type)),
             fill = "feature no.")
    }
    
    if(model == "GP") {
      # Fit
      result_filename <- file.path(output_dir,
                                   paste0("results_", use_method, "_", use_result_type, ".rds"))
      if(file.exists(result_filename)) {
        res_obj <- readRDS(result_filename)
        test_response <- res_obj$test_response
        prediction <- res_obj$prediction
        test_features <- res_obj$test_features
      } else {
        start <- Sys.time()
        gp <- mlegp(train_features, train_response)
        cat(paste0("Elapsed fit time: ", Sys.time() - start, "\n"))
  
        # Predict on test set
        start <- Sys.time()
        output_pred <- predict(gp, newData = test_features, se.fit = FALSE)
        cat(paste0("Elapsed predict time: ", Sys.time() - start, "\n"))
    
        cat("Saving results...\n")
        prediction <- output_pred[,1]
        saveRDS(list(train_features = cbind(UUID = train_uuids, train_features),
                     test_features = cbind(UUID = test_uuids, test_features),
                     train_response = train_response,
                     test_response = test_response,
                     prediction = prediction),
                result_filename)
      }
      
      plot_df <- data.frame(true = test_response,
                            predicted = prediction,
                            p = factor(test_features$P))
      levels(plot_df$p) <- c("low", "med", "high")
      pl <- ggplot(plot_df, aes(x = true, y = predicted, fill = p)) +
        geom_point(shape = 21, size = 3) +
        labs(x = paste0("observed ", toupper(use_result_type)),
             y = paste0("predicted ", toupper(use_result_type)),
             fill = "feature no.")
    }
    
    show(pl)
    ggsave(file.path("output",
                     "images",
                     paste0(model, "_results"),
                     paste0(model,
                            "_predictions_",
                            use_method,
                            "_",
                            use_result_type,
                            ".png")),
           plot = pl,
           dpi = 100,
           units = "in",
           height = 6,
           width = 7)

    # RMSE
    sq.err <- (plot_df$true - plot_df$predicted)^2
    rms.err <- sqrt(mean(sq.err))
    cat(paste0("RMSE (", model, "): ", round(rms.err, 3), "\n"))
    mean_RMSE <- mean_RMSE + rms.err
    
  }
}

cat(paste0("Mean RMSE (", model, "): ", round(mean_RMSE/10, 3), "\n"))

dbDisconnect(conn)
