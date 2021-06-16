#' Generate simulated differential expression for two conditions
#'
#' @param counts_A counts associated with the baseline/control condition 
#' (arranged p x n)
#' @param counts_B counts associated with the second condition 
#' (arranged p x n)
#' @import tidyverse
#' @import entropy
#' @import driver
#' @import moments
#' @export
characterize_dataset <- function(counts_A, counts_B) {
  if(ncol(counts_A) != ncol(counts_B)) {
    stop("Different number of features in condition A vs. B!")
  }
  p <- ncol(counts_A)
  n_a <- nrow(counts_A)
  n_b <- nrow(counts_B)
  
  pseudocount <- 0.1
  
  # Eliminate all-zero features (for now) by spiking in ones
  if(min(colSums(counts_A)) == 0) {
    counts_A <- spike_in_ones(counts_A)
  }
  if(min(colSums(counts_B)) == 0) {
    counts_B <- spike_in_ones(counts_B)
  }
  
  relab_A <- t(apply(counts_A, 1, function(x) x/sum(x)))
  relab_B <- t(apply(counts_B, 1, function(x) x/sum(x)))
  clr_A <- clr_array(counts_A + pseudocount, parts = 2)
  clr_B <- clr_array(counts_B + pseudocount, parts = 2)
  
  # --------------------------------------------------------------------------
  #   Characteristics assoc. with totals
  # --------------------------------------------------------------------------
  
  totals_A <- rowSums(counts_A)
  totals_B <- rowSums(counts_B)
  all_totals <- c(totals_A, totals_B)
  
  mean_total_A <- mean(totals_A)
  mean_total_B <- mean(totals_B)
  
  TOTALS_C_FC <- max(mean_total_B, mean_total_A) / min(mean_total_B, mean_total_A)
  TOTALS_C_D <- mean_total_B - mean_total_A
  TOTALS_C_MAX_D <- max(all_totals)
  TOTALS_C_MED_D <- median(all_totals)
  TOTALS_C_SD_D <- sd(all_totals)
  
  # --------------------------------------------------------------------------
  #   Characteristics assoc. with correlation
  # --------------------------------------------------------------------------
  
  ra_correlation <- cor(relab_A)
  ra_correlation_vector <- ra_correlation[upper.tri(ra_correlation,
                                                    diag = FALSE)]
  
  log_correlation <- cor(log(counts_A + pseudocount))
  log_correlation_vector <- log_correlation[upper.tri(log_correlation,
                                                      diag = FALSE)]
  clr_correlation <- cor(clr_A)
  clr_correlation_vector <- clr_correlation[upper.tri(clr_correlation,
                                                      diag = FALSE)]
  CORR_RA_MED <- median(ra_correlation_vector)
  CORR_RA_SD <- sd(ra_correlation_vector)
  CORR_RA_SKEW <- skewness(ra_correlation_vector)
  
  CORR_LOG_MED <- median(log_correlation_vector)
  CORR_LOG_SD <- sd(log_correlation_vector)
  CORR_LOG_SKEW <- skewness(log_correlation_vector)
  
  CORR_CLR_MED <- median(clr_correlation_vector)
  CORR_CLR_SD <- sd(clr_correlation_vector)
  CORR_CLR_SKEW <- skewness(clr_correlation_vector)
  
  # --------------------------------------------------------------------------
  #   Characteristics assoc. with composition
  # --------------------------------------------------------------------------
  
  COMP_C_P0_A <- sum(counts_A == 0) / (n_a*p)
  COMP_C_P0_B <- sum(counts_B == 0) / (n_b*p)
  COMP_C_P1_A <- sum(counts_A == 1) / (n_a*p)
  COMP_C_P1_B <- sum(counts_B == 1) / (n_b*p)
  COMP_C_P5_A <- sum(counts_A <= 5) / (n_a*p)
  COMP_C_P5_B <- sum(counts_B <= 5) / (n_b*p)
  
  COMP_RA_P01_A <- sum(relab_A < 0.001) / (n_a*p)
  COMP_RA_P01_B <- sum(relab_B < 0.001) / (n_b*p)
  COMP_RA_P1_A <- sum(relab_A < 0.01) / (n_a*p)
  COMP_RA_P1_B <- sum(relab_B < 0.01) / (n_b*p)
  COMP_RA_P5_A <- sum(relab_A < 0.05) / (n_a*p)
  COMP_RA_P5_B <- sum(relab_B < 0.05) / (n_b*p)
  
  COMP_RA_MAX_A <- max(relab_A)
  COMP_RA_MED_A <- median(relab_A)
  COMP_RA_SD_A <- sd(relab_A)
  COMP_RA_SKEW_A <- skewness(c(relab_A))
  
  COMP_RA_MAX_B <- max(relab_B)
  COMP_RA_MED_B <- median(relab_B)
  COMP_RA_SD_B <- sd(relab_B)
  COMP_RA_SKEW_B <- skewness(c(relab_B))
  
  bins <- quantile(c(relab_A[relab_A > 0], relab_B[relab_B > 0]),
                   probs = seq(from = 0, to = 1, length.out = 20))
  
  COMP_C_ENT_A <- entropy(table(cut(relab_A, breaks = c(0, bins))))
  COMP_C_ENT_B <- entropy(table(cut(relab_B, breaks = c(0, bins))))
  
  # --------------------------------------------------------------------------
  #   Characteristics assoc. with feature-wise change
  # --------------------------------------------------------------------------
  
  mean_relab_A <- colMeans(relab_A)
  mean_relab_B <- colMeans(relab_B)
  
  relab_delta <- mean_relab_B - mean_relab_A
  relab_fc <- mean_relab_B / mean_relab_A
  
  FW_RA_MAX_D <- max(relab_delta)
  FW_RA_MED_D <- median(relab_delta)
  FW_RA_SD_D <- sd(relab_delta)
  FW_RA_PPOS_D <- sum(relab_delta > 0) / p
  FW_RA_PNEG_D <- sum(relab_delta < 0) / p
  FW_RA_PFC05_D <- sum(relab_fc < 0.5) / p
  FW_RA_PFC1_D <- sum(relab_fc < 1) / p
  FW_RA_PFC2_D <- sum(relab_fc < 2) / p
  
  mean_log_A <- colMeans(log(counts_A + pseudocount))
  mean_log_B <- colMeans(log(counts_B + pseudocount))
  
  log_delta <- mean_log_B - mean_log_A
  log_fc <- mean_log_B / mean_log_A
  
  FW_LOG_MAX_D <- max(log_delta)
  FW_LOG_MED_D <- median(log_delta)
  FW_LOG_SD_D <- sd(log_delta)
  FW_LOG_PPOS_D <- sum(log_delta > 0) / p
  FW_LOG_PNEG_D <- sum(log_delta < 0) / p
  FW_LOG_PFC05_D <- sum(log_fc < 0.5) / p
  FW_LOG_PFC1_D <- sum(log_fc < 1) / p
  FW_LOG_PFC2_D <- sum(log_fc < 2) / p
  
  mean_clr_A <- colMeans(clr_A)
  mean_clr_B <- colMeans(clr_B)
  
  clr_delta <- mean_clr_B - mean_clr_A
  clr_fc <- mean_clr_B / mean_clr_A
  
  FW_CLR_MAX_D <- max(clr_delta)
  FW_CLR_MED_D <- median(clr_delta)
  FW_CLR_SD_D <- sd(clr_delta)
  FW_CLR_PPOS_D <- sum(clr_delta > 0) / p
  FW_CLR_PNEG_D <- sum(clr_delta < 0) / p
  FW_CLR_PFC05_D <- sum(clr_fc < 0.5) / p
  FW_CLR_PFC1_D <- sum(clr_fc < 1) / p
  FW_CLR_PFC2_D <- sum(clr_fc < 2) / p
  
  return(list("TOTALS_C_FC" = TOTALS_C_FC,
              "TOTALS_C_D" = TOTALS_C_D,
              "TOTALS_C_MAX_D" = TOTALS_C_MAX_D,
              "TOTALS_C_MED_D" = TOTALS_C_MED_D,
              "TOTALS_C_SD_D" = TOTALS_C_SD_D,
              "CORR_RA_MED" = CORR_RA_MED,
              "CORR_RA_SD" = CORR_RA_SD,
              "CORR_RA_SKEW" = CORR_RA_SKEW,
              "CORR_LOG_MED" = CORR_LOG_MED,
              "CORR_LOG_SD" = CORR_LOG_SD,
              "CORR_LOG_SKEW" = CORR_LOG_SKEW,
              "CORR_CLR_MED" = CORR_CLR_MED,
              "CORR_CLR_SD" = CORR_CLR_SD,
              "CORR_CLR_SKEW" = CORR_CLR_SKEW,
              "COMP_C_P0_A" = COMP_C_P0_A,
              "COMP_C_P0_B" = COMP_C_P0_B,
              "COMP_C_P1_A" = COMP_C_P1_A,
              "COMP_C_P1_B" = COMP_C_P1_B,
              "COMP_C_P5_A" = COMP_C_P5_A,
              "COMP_C_P5_B" = COMP_C_P5_B,
              "COMP_RA_P01_A" = COMP_RA_P01_A,
              "COMP_RA_P01_B" = COMP_RA_P01_B,
              "COMP_RA_P1_A" = COMP_RA_P1_A,
              "COMP_RA_P1_B" = COMP_RA_P1_B,
              "COMP_RA_P5_A" = COMP_RA_P5_A,
              "COMP_RA_P5_B" = COMP_RA_P5_B,
              "COMP_RA_MAX_A" = COMP_RA_MAX_A,
              "COMP_RA_MED_A" = COMP_RA_MED_A,
              "COMP_RA_SD_A" = COMP_RA_SD_A,
              "COMP_RA_SKEW_A" = COMP_RA_SKEW_A,
              "COMP_RA_MAX_B" = COMP_RA_MAX_B,
              "COMP_RA_MED_B" = COMP_RA_MED_B,
              "COMP_RA_SD_B" = COMP_RA_SD_B,
              "COMP_RA_SKEW_B" = COMP_RA_SKEW_B,
              "COMP_C_ENT_A" = COMP_C_ENT_A,
              "COMP_C_ENT_B" = COMP_C_ENT_B,
              "FW_RA_MAX_D" = FW_RA_MAX_D,
              "FW_RA_MED_D" = FW_RA_MED_D,
              "FW_RA_SD_D" = FW_RA_SD_D,
              "FW_RA_PPOS_D" = FW_RA_PPOS_D,
              "FW_RA_PNEG_D" = FW_RA_PNEG_D,
              "FW_RA_PFC05_D" = FW_RA_PFC05_D,
              "FW_RA_PFC1_D" = FW_RA_PFC1_D,
              "FW_RA_PFC2_D" = FW_RA_PFC2_D,
              "FW_LOG_MAX_D" = FW_LOG_MAX_D,
              "FW_LOG_MED_D" = FW_LOG_MED_D,
              "FW_LOG_SD_D" = FW_LOG_SD_D,
              "FW_LOG_PPOS_D" = FW_LOG_PPOS_D,
              "FW_LOG_PNEG_D" = FW_LOG_PNEG_D,
              "FW_LOG_PFC05_D" = FW_LOG_PFC05_D,
              "FW_LOG_PFC1_D" = FW_LOG_PFC1_D,
              "FW_LOG_PFC2_D" = FW_LOG_PFC2_D,
              "FW_CLR_MAX_D" = FW_CLR_MAX_D,
              "FW_CLR_MED_D" = FW_CLR_MED_D,
              "FW_CLR_SD_D" = FW_CLR_SD_D,
              "FW_CLR_PPOS_D" = FW_CLR_PPOS_D,
              "FW_CLR_PNEG_D" = FW_CLR_PNEG_D,
              "FW_CLR_PFC05_D" = FW_CLR_PFC05_D,
              "FW_CLR_PFC1_D" = FW_CLR_PFC1_D,
              "FW_CLR_PFC2_D" = FW_CLR_PFC2_D))
}
  
#' Generate simulated differential expression for two conditions
#'
#' @param model options are "linear", "RF" (random forest), "EN" (elastic net),
#' "GP" (Gaussian process)
#' @param do_predict flag indicating whether or not to generated predictions on a
#' test set
#' @param DE_method which DE calling method's results to predict; if "all",
#' prediction is over all results together
#' @param plot_weights flag indicating whether or not to plot some visualization
#' of feature weights
#' @return fitted model
#' @import RSQLite
#' @import mlegp
#' @import randomForest
#' @import tidyr
#' @import caret
#' @import jtools
#' @import ggstance
#' @import broom.mixed
#' @export
fit_predictive_model <- function(model = "linear", do_predict = FALSE,
                                 DE_method = "all", plot_weights = FALSE) {
  if(!(model %in% c("linear", "RF", "EN", "GP"))) {
    stop(paste0("Invalid model type: ", model, "!"))
  }
  if(!(DE_method %in% c("all", "ALDEx2", "DESeq2", "MAST", "NBGLM", "scran"))) {
    stop(paste0("Invalid DE calling method: ", DE_method, "!\n"))
  }
  # GP models take a while to run. Save the results for later tinkering.
  output_dir <- "output"
  if(model == "GP") {
    if(model == "GP") {
      output_dir <- file.path("output", model)
      if(!dir.exists(output_dir)) {
        dir.create(output_dir)
      }
    }
  }
  if(plot_weights) {
    result_dir <- file.path(output_dir, "images", paste0(model, "_results"))
    if(!dir.exists(result_dir)) {
      dir.create(result_dir)
    }
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
  dbDisconnect(conn)
  
  fitted_models <- list()
  predictions <- list()

  for(use_result_type in c("fpr", "tpr")) {
    cat(paste0("Modeling ", use_result_type, " w/ DE method ", DE_method, "\n"))
    
    if(DE_method == "all") {
      data <- results %>%
        filter(RESULT_TYPE == use_result_type)
    } else {
      data <- results %>%
        filter(METHOD == DE_method) %>%
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
    if(DE_method == "all") {
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
    
    factors <- which(colnames(features) %in% c("PARTIAL", "CORRP", "METHOD", "P"))
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
      res <- train(RESULT ~ ., data = train_data,
                   method = "glmnet",
                   trControl = trainControl("cv", number = 10),
                   tuneLength = 10)
      diff <- Sys.time() - start
      cat(paste0("EN model train time: ", round(diff, 2), " ", attr(diff, "units"), "\n"))
      
      # plot(res) # best parameter
      
      res$bestTune # best coefficient
      # Per the documentation:
      #   lambda seems to be the strength of the overall penalty term
      #   (1-alpha) is the ridge weight, alpha is the LASSO weight
      
      model_features <- as.matrix(coef(res$finalModel, res$bestTune$lambda))
      # Eliminate zero-weight features
      model_features <- data.frame(feature_name = rownames(model_features),
                                   beta = unname(unlist(model_features)))
      model_features <- model_features %>%
        filter(beta > 0) %>%
        filter(feature_name != "(Intercept)")
      
      if(plot_weights) {
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
                                DE_method,
                                "_",
                                use_result_type,
                                ".png")),
               plot = pl,
               dpi = 100,
               units = "in",
               height = 8,
               width = 4)
      }
      
      if(do_predict) {
        test_data <- cbind(RESULT = response_vector[-trainRowNumbers],
                           features[-trainRowNumbers,])
        
        test_response <- test_data$RESULT
        prediction <- predict(res, test_data)
        p_labels <- test_data$P
      }
    } else {
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
      
      if(model == "linear") {
        lm_train_data <- cbind(train_features, train_response)
        res <- lm(train_response ~ ., data = lm_train_data)
        if(plot_weights) {
          pl <- plot_summs(res)
          ggsave(file.path("output",
                           "images",
                           paste0(model, "_results"),
                           paste0(model,
                                  "_betas_",
                                  DE_method,
                                  "_",
                                  use_result_type,
                                  ".png")),
                 plot = pl,
                 dpi = 100,
                 units = "in",
                 height = 8,
                 width = 4)
        }
        if(do_predict) {
          lm_test_data <- cbind(test_features, test_response)
          prediction <- predict(res, newdata = lm_test_data)
          test_response <- test_response
          p_labels <- test_features$P
        }
      }
      
      if(model == "RF") {
        rf_train_data <- cbind(train_features, train_response)
        res <- randomForest(train_response ~ ., data = rf_train_data)
        if(do_predict) {
          rf_test_data <- cbind(test_features, test_response)
          prediction <- predict(res, newdata = rf_test_data)
          test_response <- test_response
          p_labels <- test_features$P
        }
      }
      
      if(model == "GP") {
        # Fit
        result_filename <- file.path(output_dir,
                                     paste0("results_", DE_method, "_", use_result_type, ".rds"))
        if(do_predict & file.exists(result_filename)) {
          res_obj <- readRDS(result_filename)
          test_response <- res_obj$test_response
          prediction <- res_obj$prediction
          test_features <- res_obj$test_features
          p_labels <- test_features$P
        } else {
          start <- Sys.time()
          res <- mlegp(train_features, train_response)
          cat(paste0("Elapsed fit time: ", Sys.time() - start, "\n"))
          
          if(do_predict) {
            # Predict on test set
            start <- Sys.time()
            output_pred <- predict(res, newData = test_features, se.fit = FALSE)
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
    }
    fitted_models[[use_result_type]] <- res
    if(do_predict) {
      predictions[[use_result_type]] <- list(true = test_response,
                                             predicted = prediction,
                                             p_labels = p_labels)
    }
  }
  if(do_predict) {
    return(list(fitted_model = fitted_models, predictions = predictions))
  } else {
    return(fitted_models)
  }
}

