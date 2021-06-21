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
                                   "METHOD, RESULT, RESULT_TYPE FROM ",
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

    # These predictors appear to be strongly correlated.
    # Let's drop them for now.
    data <- data %>%
      select(-c(FW_RA_PFC1_D, FW_CLR_MED_D, FW_CLR_SD_D, FW_CLR_PNEG_D))
    
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
    
    # Could add elastic net here!
    
    if(model == "linear") {
      lm_train_data <- cbind(train_features, train_response)
      res <- lm(train_response ~ ., data = lm_train_data)
      pl <- plot_summs(res)
      ggsave(file.path("output",
                       "images",
                       paste0(model, "_results"),
                       paste0(model,
                              "_betas_",
                              use_method,
                              "_",
                              use_result_type,
                              ".png")),
             plot = pl,
             dpi = 100,
             units = "in",
             height = 8,
             width = 4)
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
      res <- randomForest(train_response ~ ., data = rf_train_data)
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
    
    # show(pl)
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
