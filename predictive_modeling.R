source("path_fix.R")

library(codaDE)
library(tidyverse)
library(RSQLite)
library(mlegp)
library(randomForest)
library(tidyr)
library(caret)
library(optparse)

# For plotting betas
library(jtools)
library(ggstance)
library(broom.mixed)

option_list = list(
  make_option(c("--model"), type = "character", default = "EN",
              help = "predictive model to use: EN, GP, linear, RF", metavar = "character"),
  make_option(c("--together"), type = "logical", default = "FALSE",
              help = "combine results of all methods", metavar = "logical")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);
model <- opt$model
if(!(model %in% c("EN", "GP", "linear", "RF"))) {
  stop("Invalid model specified!\n")
}
together <- opt$together

# model <- "EN"
# together <- TRUE

# GP models take a while to run. Save the results for later tinkering.
output_dir <- "output"
if(model == "GP") {
  output_dir <- file.path("output", model)
  if(!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
}
result_dir <- file.path(output_dir, "images", paste0(model, "_results"))
if(!dir.exists(result_dir)) {
  dir.create(result_dir)
}

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

results <- dbGetQuery(conn, paste0("SELECT datasets.UUID, P, PARTIAL, CORRP, ",
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

# mean_RMSE <- 0

method_list <- c("ALDEx2", "DESeq2", "MAST", "NBGLM", "scran")
if(together) {
  method_list <- c("all")
}

plot_labels <- list(fpr = "specificity (1 - FPR)", tpr = "sensitivity (TPR)")

for(use_method in method_list) {
  for(use_result_type in c("fpr", "tpr")) {

    cat(paste0("Modeling ", use_result_type, " w/ ", use_method, "\n"))

    if(use_method == "all") {
      data <- results %>%
        filter(RESULT_TYPE == use_result_type)
    } else {
      data <- results %>%
        filter(METHOD == use_method) %>%
        filter(RESULT_TYPE == use_result_type)
    }
    
    # These predictors appear to be strongly correlated.
    # Let's drop them for now.
    data <- data %>%
      select(-c(FW_RA_PFC1_D, FW_CLR_MED_D, FW_CLR_SD_D, FW_CLR_PNEG_D))
    
    # Subset for testing
    # data <- data[sample(1:nrow(data), size = 100),]

    # Scale the features
    uuids <- data$UUID
    if(use_method == "all") {
      features <- data %>%
        select(!c(UUID, RESULT, RESULT_TYPE))
      features$PARTIAL <- factor(features$PARTIAL)
      features$CORRP <- factor(features$CORRP)
      features$METHOD <- factor(features$METHOD)
    } else {
      features <- data %>%
        select(!c(UUID, METHOD, RESULT, RESULT_TYPE))
      features$PARTIAL <- factor(features$PARTIAL)
      features$CORRP <- factor(features$CORRP)
    }

    factors <- which(colnames(features) %in% c("PARTIAL", "CORRP", "METHOD"))
    non_factors <- setdiff(1:ncol(features), factors)
    features_nonfactors <- features[,non_factors]
    features_nonfactors <- as.data.frame(apply(features_nonfactors, 2, function(x) {
      if(sd(x) > 0) {
        scale(x)
      } else {
        scale(x, scale = FALSE)
      }
    }))
    # Drop columns with no variation
    # In practice this only happens in testing/subsetting to small samples
    features_nonfactors <- features_nonfactors[,apply(features_nonfactors, 2, sd) > 0]
    features <- cbind(features[,factors], features_nonfactors)
    response <- data %>%
      select(RESULT)

    n <- nrow(data)

    if(model == "EN") {
      response_vector <- unname(unlist(response))
      if(use_result_type == "fpr") {
        # Predict specificity: 1 - fpr
        response_vector <- 1 - response_vector
      }
      trainRowNumbers <- createDataPartition(response_vector,
                                             p = 0.5,
                                             list = FALSE)
      
      # Separate out the training set
      train_data <- cbind(RESULT = response_vector[trainRowNumbers],
                          features[trainRowNumbers,])

      start <- Sys.time()
      model_res <- train(RESULT ~ ., data = train_data,
                     method = "glmnet",
                     trControl = trainControl("cv", number = 10),
                     tuneLength = 10)
      diff <- Sys.time() - start
      cat(paste0("EN model train time: ", round(diff, 2), " ", attr(diff, "units"), "\n"))
      
      plot(model_res) # best parameter
      
      model_res$bestTune # best coefficient
      # Per the documentation:
      #   lambda seems to be the strength of the overall penalty term
      #   (1-alpha) is the ridge weight, alpha is the LASSO weight
      
      model_features <- as.matrix(coef(model_res$finalModel, model_res$bestTune$lambda))
      # Eliminate zero-weight features
      model_features <- data.frame(feature_name = rownames(model_features),
                                   beta = unname(unlist(model_features)))
      model_features <- model_features %>%
        filter(beta > 0) %>%
        filter(feature_name != "(Intercept)")
      
      pl <- ggplot(model_features, aes(reorder(feature_name, beta), beta)) +
        geom_bar(stat = "identity") + 
        coord_flip() + 
        scale_y_continuous("Weight") +
        scale_x_discrete("Ordered feature weights")
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
      
      test_data <- cbind(RESULT = response_vector[-trainRowNumbers],
                         features[-trainRowNumbers,])
      
      test_response <- test_data$RESULT
      prediction <- predict(model_res, test_data)
      p_labels <- test_data$P
    } else {
      # Define test/train set
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
        prediction <- predict(res, newdata = lm_test_data)
        test_response <- test_response
        p_labels <- test_features$P
      }
  
      if(model == "RF") {
        rf_train_data <- cbind(train_features, train_response)
        res <- randomForest(train_response ~ ., data = rf_train_data)
        rf_test_data <- cbind(test_features, test_response)
        prediction <- predict(res, newdata = rf_test_data)
        test_response <- test_response
        p_labels <- test_features$P
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
      }
    }

    plot_df <- data.frame(true = test_response,
                          predicted = prediction,
                          p = factor(p_labels))
    levels(plot_df$p) <- c("low", "med", "high")
    pl <- ggplot(plot_df, aes(x = true, y = predicted, fill = p)) +
      geom_point(shape = 21, size = 3) +
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      labs(x = paste0("observed ", plot_labels[[use_result_type]]),
           y = paste0("predicted ", plot_labels[[use_result_type]]),
           fill = "feature no.")

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
           height = 4,
           width = 5)

    # RMSE
    # sq.err <- (plot_df$true - plot_df$predicted)^2
    # rms.err <- sqrt(mean(sq.err))
    # cat(paste0("RMSE (", model, "): ", round(rms.err, 3), "\n"))
    # mean_RMSE <- mean_RMSE + rms.err
    
    cat(paste0("R^2:", round(cor(test_response, prediction)^2, 3), "\n"))
  }
}

# cat(paste0("Mean RMSE (", model, "): ", round(mean_RMSE/10, 3), "\n"))

dbDisconnect(conn)
