#' Generate simulated differential expression for two conditions
#'
#' @param counts_A counts associated with the baseline/control condition 
#' (arranged n x p)
#' @param counts_B counts associated with the second condition 
#' (arranged n x p)
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
  
  temp <- log(counts_A + pseudocount)
  # Find features with zero standard deviation and remove them here
  # They'll cause trouble with the correlation computation
  remove_idx <- unname(which(apply(temp, 2, sd) == 0))
  if(length(remove_idx) > 0) {
    temp <- temp[,-remove_idx]
  }
  log_correlation <- cor(temp)
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
  # Add jitter to force unique bins
  bins <- bins + seq_along(bins) * .Machine$double.eps

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
#' @param DE_method which DE calling method's results to predict; if "all",
#' prediction is over all results together
#' @param use_baseline "self" or "oracle"
#' @param plot_weights flag indicating whether or not to plot some visualization
#' of feature weights
#' @param exclude_partials flag indicating whether to exclude simulations with 
#' partial abundance information
#' @param exclude_independent flag indicating whether to exclude simulations 
#' with uncorrelated features
#' @param train_percent percent of simulated datasets to train on
#' @param save_slug optional file path and name for saved predictive model (e.g.
#' /output/allmethods_oracle, to which _TPR.rds and _FPR.rds will be appended
#' @param save_training_data flag indicating whether or not to save training data
#' with the saved model output
#' @return NULL (fitted models are saved in output directory)
#' @import RSQLite
#' @import randomForest
#' @import tidyr
#' @import caret
#' @import jtools
#' @import ggstance
#' @import broom.mixed
#' @export
fit_predictive_model <- function(DE_method = "all",
                                 use_baseline = "self",
                                 plot_weights = FALSE,
                                 exclude_partials = TRUE,
                                 exclude_independent = FALSE,
                                 train_percent = 0.8,
                                 save_slug = NULL,
                                 save_training_data = TRUE) {
  if(!(use_baseline %in% c("self", "oracle"))) {
    stop(paste0("Invalid baseline: ", use_baseline, "!\n"))
  }

  if(!(DE_method %in% c("all", "ALDEx2", "DESeq2", "MAST", "scran"))) {
    stop(paste0("Invalid DE calling method: ", DE_method, "!\n"))
  }

  save_path <- list("output", "predictive_fits", DE_method)
  for(i in 1:length(save_path)) {
    fp <- do.call(file.path, save_path[1:i])
    if(!dir.exists(fp)) {
      dir.create(fp)
    }
  }
  save_dir <- fp
  
  conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

  results <- dbGetQuery(conn, paste0("SELECT datasets.UUID AS UUID, P, CORRP, ",
                                     "FC_ABSOLUTE, FC_PARTIAL, FC_RELATIVE, ",
                                     "TOTALS_C_FC, TOTALS_C_D, ",
                                     "TOTALS_C_MAX_D, TOTALS_C_MED_D, ",
                                     "TOTALS_C_SD_D, CORR_RA_MED, CORR_RA_SD, ",
                                     "CORR_RA_SKEW, CORR_LOG_MED, CORR_LOG_SD, ",
                                     "CORR_LOG_SKEW, CORR_CLR_MED, CORR_CLR_SD, ",
                                     "CORR_CLR_SKEW, COMP_C_P0_A, COMP_C_P0_B, ",
                                     "COMP_C_P1_A, COMP_C_P1_B, COMP_C_P5_A, ",
                                     "COMP_C_P5_B, COMP_RA_P01_A, COMP_RA_P01_B, ",
                                     "COMP_RA_P1_A, COMP_RA_P1_B, COMP_RA_P5_A, ",
                                     "COMP_RA_P5_B, COMP_RA_MAX_A, COMP_RA_MED_A, ",
                                     "COMP_RA_SD_A, COMP_RA_SKEW_A, COMP_RA_MAX_B, ",
                                     "COMP_RA_MED_B, COMP_RA_SD_B, COMP_RA_SKEW_B, ",
                                     "COMP_C_ENT_A, COMP_C_ENT_B, FW_RA_MAX_D, ",
                                     "FW_RA_MED_D, FW_RA_SD_D, FW_RA_PPOS_D, ",
                                     "FW_RA_PNEG_D, FW_RA_PFC05_D, FW_RA_PFC1_D, ",
                                     "FW_RA_PFC2_D, FW_LOG_MAX_D, FW_LOG_MED_D, ",
                                     "FW_LOG_SD_D, FW_LOG_PPOS_D, FW_LOG_PNEG_D, ",
                                     "FW_LOG_PFC05_D, FW_LOG_PFC1_D, ",
                                     "FW_LOG_PFC2_D, FW_CLR_MAX_D, FW_CLR_MED_D, ",
                                     "FW_CLR_SD_D, FW_CLR_PPOS_D, FW_CLR_PNEG_D, ",
                                     "FW_CLR_PFC05_D, FW_CLR_PFC1_D, ",
                                     "FW_CLR_PFC2_D, METHOD, PARTIAL_INFO, ",
                                     "BASELINE_TYPE,TPR, FPR ",
                                     "FROM datasets LEFT JOIN characteristics ",
                                     "ON datasets.UUID=characteristics.UUID ",
                                     "LEFT JOIN results ",
                                     "ON (characteristics.UUID=results.UUID ",
                                     "AND characteristics.partial=results.PARTIAL_INFO);"))
  
  dbDisconnect(conn)
  
  if(exclude_partials) {
    # Filter out results with partial information
    results <- results %>%
      filter(PARTIAL_INFO == 0)
  } else {
    # Filter simulations with partial abundances to those that are "informative"
    # i.e. those with an a fold change intermediate between the observed FC (~1)
    # and the true FC
    results$include <- FALSE
    results$include[results$PARTIAL_INFO == 0] <- TRUE
    # Decreases
    results$include[(results$PARTIAL_INFO == 1 &
                       results$FC_ABSOLUTE > 0.5 &
                       results$FC_ABSOLUTE < 1 &
                       results$FC_ABSOLUTE < results$FC_PARTIAL &
                       results$FC_PARTIAL < results$FC_RELATIVE)] <- TRUE
    # Increases
    results$include[(results$PARTIAL_INFO == 1 &
                       results$FC_ABSOLUTE < 2 &
                       results$FC_ABSOLUTE > 1 &
                       results$FC_ABSOLUTE > results$FC_PARTIAL &
                       results$FC_PARTIAL > results$FC_RELATIVE)] <- TRUE
    results <- results %>%
      filter(include) %>%
      select(!c("include"))
  }
  
  results <- results %>%
    select(!c("FC_ABSOLUTE", "FC_PARTIAL", "FC_RELATIVE", "PARTIAL_INFO"))
  
  if(exclude_independent) {
    # Filter out uncorrelated-feature simulations
    results <- results %>%
      filter(CORRP != 0)
  }
  
  # Strip "result-less" entries
  results <- results %>%
    filter(!is.na(TPR) & !is.na(FPR))
  
  for(use_result_type in c("TPR", "FPR")) {
    cat(paste0("Modeling ", use_result_type, " w/ DE method ", DE_method, "\n"))
    
    data <- results %>%
      filter(BASELINE_TYPE == use_baseline)
    
    # These predictors appear to be strongly correlated.
    # Let's drop them for now.
    data <- data %>%
      select(-c(FW_RA_PFC1_D, FW_CLR_MED_D, FW_CLR_SD_D, FW_CLR_PNEG_D))
    
    # Pull out the features we won't predict on
    uuids <- data$UUID
    if(DE_method == "all") {
      features <- data %>%
        select(-c(UUID, CORRP, BASELINE_TYPE))
      features$METHOD <- factor(features$METHOD)
    } else {
      features <- data %>%
        filter(METHOD == DE_method) %>%
        select(-c(UUID, CORRP, BASELINE_TYPE, METHOD))
    }
    
    # Remove result type we're not interested in and separate data into a
    # features and response data.frame
    features <- features %>%
      select(-one_of(ifelse(use_result_type == "TPR", "FPR", "TPR")))
    response <- features %>%
      select(one_of(use_result_type))
    features <- features %>%
      select(-one_of(use_result_type))
    
    # factors <- which(colnames(features) %in% c("METHOD", "P"))
    # non_factors <- setdiff(1:ncol(features), factors)
    # features_nonfactors <- features[,non_factors]
    # Drop columns with no variation
    # In practice this only happens in testing/subsetting to small samples
    # features_nonfactors <- features_nonfactors[,apply(features_nonfactors, 2, sd) > 0]
    # features <- cbind(features[,factors,drop=F], features_nonfactors)
    
    n <- nrow(features)
    
    # Define test/train set for all remaining methods
    train_idx <- sample(1:n, size = round(n*train_percent))
    test_idx <- setdiff(1:n, train_idx)
    train_uuids <- uuids[train_idx]
    train_features <- features[train_idx,]
    train_response <- unname(unlist(response[train_idx,1]))
    if(use_result_type == "FPR") {
      train_response <- 1 - train_response
    }
    test_uuids <- uuids[test_idx]
    test_features <- features[test_idx,]
    test_response <- unname(unlist(response[test_idx,1]))
    if(use_result_type == "FPR") {
      test_response <- 1 - test_response
    }

    if(is.null(save_slug)) {
      save_fn <- file.path(save_dir,
                           paste0(DE_method,
                                  "_",
                                  use_baseline,
                                  "_",
                                  use_result_type,
                                  ".rds"))
    } else {
      save_fn <- paste0(save_slug, "_", use_result_type, ".rds")
    }
    cat(paste0("Predictive model saving/saved to: ", save_fn, "\n"))
    if(!file.exists(save_fn)) {
      cat("Training model...\n")
      rf_train_data <- cbind(train_features, response = train_response)
      res <- randomForest(response ~ ., data = rf_train_data)
      if(save_training_data) {
        saveRDS(list(result = res,
                     train_features = train_features,
                     train_response = train_response,
                     test_features = test_features,
                     test_response = test_response), save_fn)
      } else {
        saveRDS(list(result = res), save_fn)
      }
    } else {
      cat("Loading existing model...\n")
      res_obj <- readRDS(save_fn)
      res <- res_obj$result
      train_features <- res_obj$train_features
      train_response <- res_obj$train_response
      test_features <- res_obj$test_features
      test_response <- res_obj$test_response
    }

    if(plot_weights) {
      # Need to find out if there's a similar function for rfinterval
      plot_df <- varImpPlot(res)
      plot_df <- data.frame(feature = rownames(plot_df),
                            "IncNodePurity" = plot_df[,1])
      rownames(plot_df) <- NULL
      plot_df <- plot_df %>%
        arrange(desc(IncNodePurity)) %>%
        slice(1:20)
      pl <- ggplot(plot_df, aes(x = IncNodePurity,
                                y = reorder(feature, IncNodePurity))) +
        geom_bar(stat = "identity") +
        labs(y = "feature")
      # show(pl)
      ggsave(file.path(save_dir,
                       paste0("betas_",
                              DE_method,
                              "_",
                              use_result_type,
                              "_",
                              use_baseline,
                              ".png")),
             plot = pl,
             dpi = 100,
             units = "in",
             height = 4,
             width = 6)
    }
    
    # To predict do:
    # rf_test_data <- cbind(test_features, test_response)
    # prediction <- predict(res, newdata = rf_test_data)
  }
}

