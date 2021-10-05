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
  
  if(ncol(relab_A) > 10000) {
    ra_correlation <- cor(relab_A[,1:10000])
  } else {
    ra_correlation <- cor(relab_A)
  }
  ra_correlation_vector <- ra_correlation[upper.tri(ra_correlation,
                                                    diag = FALSE)]
  
  temp <- log(counts_A + pseudocount)
  # Find features with zero standard deviation and remove them here
  # They'll cause trouble with the correlation computation
  remove_idx <- unname(which(apply(temp, 2, sd) == 0))
  if(length(remove_idx) > 0) {
    temp <- temp[,-remove_idx]
  }
  if(ncol(temp) > 10000) {
    log_correlation <- cor(temp[,1:10000])
  } else {
    log_correlation <- cor(temp)
  }
  log_correlation_vector <- log_correlation[upper.tri(log_correlation,
                                                      diag = FALSE)]

  if(ncol(clr_A) > 10000) {
    clr_correlation <- cor(clr_A[,1:10000])
  } else {
    clr_correlation <- cor(clr_A)
  }
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
#' @param DE_methods differential abundance calling methods to use
#' @param use_baseline one of either "self" or "oracle"; baseline differential
#' abundance calls against which we will score accuracy of calls made on 
#' observed (relative) abundance data
#' @param exclude_partials flag indicating whether to exclude simulations with 
#' partial abundance information
#' @param exclude_independent flag indicating whether to exclude simulations 
#' with uncorrelated features
#' @param feature_list optional predictive feature list to use
#' @param abs_feature_list optional predictive feature list from absolute 
#' abundance data to use
#' @return data.frame with predictive model training features
#' @import dplyr
#' @import RSQLite
#' @export
pull_features <- function(DE_methods = c("ALDEx2", "DESeq2", "scran"),
                          use_baseline = "self",
                          exclude_partials = TRUE,
                          exclude_independent = FALSE,
                          feature_list = NULL,
                          abs_feature_list = NULL) {

  conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
  
  # Note: I'm not pulling FW_RA_PFC1_D, FW_CLR_MED_D, FW_CLR_SD_D, FW_CLR_PNEG_D
  # because these predictors appear to be strongly correlated.
  
  if(!is.null(feature_list) & length(feature_list > 0)) {
    if(!is.null(abs_feature_list) & length(abs_feature_list > 0)) {
      results <- dbGetQuery(conn, paste0("SELECT ", paste(unique(c(feature_list, abs_feature_list)), collapse = ", "), " ",
                                         "FROM datasets LEFT JOIN characteristics ",
                                         "ON datasets.UUID=characteristics.UUID ",
                                         "LEFT JOIN results ",
                                         "ON (characteristics.UUID=results.UUID ",
                                         "AND characteristics.partial=results.PARTIAL_INFO) ",
                                         "WHERE BASELINE_TYPE='", use_baseline, "';"))
    } else {
      results <- dbGetQuery(conn, paste0("SELECT ", paste(feature_list, collapse = ", "), " ",
                                       "FROM datasets LEFT JOIN characteristics ",
                                       "ON datasets.UUID=characteristics.UUID ",
                                       "LEFT JOIN results ",
                                       "ON (characteristics.UUID=results.UUID ",
                                       "AND characteristics.partial=results.PARTIAL_INFO) ",
                                       "WHERE BASELINE_TYPE='", use_baseline, "';"))
    }
    dbDisconnect(conn)
    return(results)
  } else {
    if(!is.null(abs_feature_list) & length(abs_feature_list > 0)) {
      results <- dbGetQuery(conn, paste0("SELECT datasets.UUID AS UUID, P, CORRP, ",
                                         paste0(paste(abs_feature_list, collapse = ", "), ", "),
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
                                         "FW_RA_PNEG_D, FW_RA_PFC05_D, ",
                                         "FW_RA_PFC2_D, FW_LOG_MAX_D, FW_LOG_MED_D, ",
                                         "FW_LOG_SD_D, FW_LOG_PPOS_D, FW_LOG_PNEG_D, ",
                                         "FW_LOG_PFC05_D, FW_LOG_PFC1_D, ",
                                         "FW_LOG_PFC2_D, FW_CLR_MAX_D, ",
                                         "FW_CLR_PPOS_D, ",
                                         "FW_CLR_PFC05_D, FW_CLR_PFC1_D, ",
                                         "FW_CLR_PFC2_D, METHOD, PARTIAL_INFO, ",
                                         "BASELINE_TYPE,TPR, FPR ",
                                         "FROM datasets LEFT JOIN characteristics ",
                                         "ON datasets.UUID=characteristics.UUID ",
                                         "LEFT JOIN results ",
                                         "ON (characteristics.UUID=results.UUID ",
                                         "AND characteristics.partial=results.PARTIAL_INFO) ",
                                         "WHERE BASELINE_TYPE='", use_baseline, "';"))
    } else {
      results <- dbGetQuery(conn, paste0("SELECT datasets.UUID AS UUID, P, CORRP, ",
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
                                         "FW_RA_PNEG_D, FW_RA_PFC05_D, ",
                                         "FW_RA_PFC2_D, FW_LOG_MAX_D, FW_LOG_MED_D, ",
                                         "FW_LOG_SD_D, FW_LOG_PPOS_D, FW_LOG_PNEG_D, ",
                                         "FW_LOG_PFC05_D, FW_LOG_PFC1_D, ",
                                         "FW_LOG_PFC2_D, FW_CLR_MAX_D, ",
                                         "FW_CLR_PPOS_D, ",
                                         "FW_CLR_PFC05_D, FW_CLR_PFC1_D, ",
                                         "FW_CLR_PFC2_D, METHOD, PARTIAL_INFO, ",
                                         "BASELINE_TYPE,TPR, FPR ",
                                         "FROM datasets LEFT JOIN characteristics ",
                                         "ON datasets.UUID=characteristics.UUID ",
                                         "LEFT JOIN results ",
                                         "ON (characteristics.UUID=results.UUID ",
                                         "AND characteristics.partial=results.PARTIAL_INFO) ",
                                         "WHERE BASELINE_TYPE='", use_baseline, "';"))
    }
    dbDisconnect(conn)
  }

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
    select(!PARTIAL_INFO)
  
  if(exclude_independent) {
    # Filter out uncorrelated-feature simulations
    results <- results %>%
      filter(CORRP != 0)
  }
  
  # Strip "result-less" entries
  results <- results %>%
    filter(!is.na(TPR) & !is.na(FPR))
  
  # Pull out the features we won't predict on
  uuids <- results$UUID
  results <- results %>%
    filter(METHOD %in% DE_methods) %>%
    select(-c(UUID, CORRP, BASELINE_TYPE))
  
  # With UUID
  # results <- results %>%
  #   filter(METHOD %in% DE_methods) %>%
  #   select(-c(CORRP, BASELINE_TYPE))
  
  return(results)
}

#' Generate simulated differential expression for two conditions
#'
#' @param DE_methods which DE calling method's results to predict; if "all",
#' prediction is over all results together
#' @param use_baseline "self" or "oracle"
#' @param output_weights flag indicating whether or not to plot some visualization
#' of feature weights
#' @param exclude_partials flag indicating whether to exclude simulations with 
#' partial abundance information
#' @param exclude_independent flag indicating whether to exclude simulations 
#' with uncorrelated features
#' @param train_percent percent of simulated datasets to train on
#' @param abs_feature_list optional predictive feature list from absolute 
#' abundance data to use
#' @param save_slug optional file path and name for saved predictive model (e.g.
#' /output/allmethods_oracle, to which _TPR.rds and _FPR.rds will be appended
#' @param save_training_data flag indicating whether or not to save training data
#' with the saved model output
#' @param do_classify optional flag indicating whether or not to train a 
#' classification model on "good" vs. "bad" outcomes
#' @param alpha optional threshold for classification indicating the minimum 
#' sensitivity of specificity for "good" outcomes
#' @return NULL (fitted models are saved in output directory)
#' @import randomForest
#' @import dplyr
#' @import tidyr
#' @import caret
#' @export
fit_predictive_model <- function(DE_methods = c("ALDEx2", "DESeq2", "scran"),
                                 use_baseline = "self",
                                 output_weights = TRUE,
                                 exclude_partials = TRUE,
                                 exclude_independent = FALSE,
                                 train_percent = 0.8,
                                 abs_feature_list = NULL,
                                 save_slug = NULL,
                                 save_training_data = TRUE,
                                 do_classify = FALSE,
                                 alpha = 0.95) {
  
  if(!(use_baseline %in% c("self", "oracle"))) {
    stop(paste0("Invalid baseline: ", use_baseline, "!\n"))
  }
  
  if(!any(DE_methods %in% c("ALDEx2", "DESeq2", "MAST", "scran"))) {
    stop(paste0("Invalid DE calling method!\n"))
  }

  save_path <- list("output",
                    "predictive_fits",
                    paste0(use_baseline, "_", ifelse(exclude_partials, "nopartial", "partial")),
                    ifelse(do_classify, "classification", "regression"))
  for(i in 1:length(save_path)) {
    fp <- do.call(file.path, save_path[1:i])
    if(!dir.exists(fp)) {
      dir.create(fp)
    }
  }
  save_dir <- fp
  
  features <- pull_features(use_baseline = use_baseline,
                            exclude_partials = exclude_partials,
                            exclude_independent = exclude_independent,
                            abs_feature_list = abs_feature_list)
  
  # Subset to methods of interest
  if(length(DE_methods) > 1) {
    features <- features %>%
      filter(METHOD %in% DE_methods)
    features$METHOD <- factor(features$METHOD)
  } else {
    features <- features %>%
      filter(METHOD %in% DE_methods) %>%
      select(!METHOD)
  }
  
  # for(use_result_type in c("TPR", "FPR")) {
  for(use_result_type in c("TPR", "FPR")) {
    cat(paste0("Modeling ", use_result_type, "\n"))
    
    # Pull out the features we won't predict on
    uuids <- features$UUID

    # Remove result type we're not interested in and separate data into a
    # features and response data.frame
    features_result_type <- features %>%
      select(-one_of(ifelse(use_result_type == "TPR", "FPR", "TPR")))
    response <- features_result_type %>%
      select(one_of(use_result_type))
    features_result_type <- features_result_type %>%
      select(-one_of(use_result_type))
    
    # Subset for testing
    # idx <- sample(1:nrow(features), size = 1000)
    # features_result_type <- features_result_type[idx,]
    # response <- response[idx,,drop=FALSE]
    
    n <- nrow(features_result_type)
    
    # Define test/train set for all remaining methods
    train_idx <- sample(1:n, size = round(n*train_percent))
    test_idx <- setdiff(1:n, train_idx)
    train_uuids <- uuids[train_idx]
    train_features <- features_result_type[train_idx,]
    train_response <- unname(unlist(response[train_idx,1]))
    if(use_result_type == "FPR") {
      train_response <- 1 - train_response
    }
    test_uuids <- uuids[test_idx]
    test_features <- features_result_type[test_idx,]
    test_response <- unname(unlist(response[test_idx,1]))
    if(use_result_type == "FPR") {
      test_response <- 1 - test_response
    }
    
    if(do_classify) {
      train_response <- sapply(train_response, function(x) {
        # "Bad" outcomes (1) have less than 95% sensitivity or specificity
        # "Good" outcomes (0) have at least 95% sensitivity or specificity
        ifelse(x < alpha, 1, 0)
      })
      train_response <- factor(train_response)
      test_response <- sapply(test_response, function(x) {
        ifelse(x < alpha, 1, 0)
      })
      test_response <- factor(test_response)
    }

    if(is.null(save_slug)) {
      save_fn <- file.path(save_dir,
                           paste0(use_result_type,
                                  ifelse(length(DE_methods) == 1, paste0("_", DE_methods), ""),
                                  ".rds"))
    } else {
      save_fn <- file.path(save_dir,
                           paste0(save_slug,
                                  "_",
                                  use_result_type,
                                  ifelse(length(DE_methods) == 1, paste0("_", DE_methods), ""),
                                  ".rds"))
    }
    cat(paste0("Predictive model saving/saved to: ", save_fn, "\n"))
    if(!file.exists(save_fn)) {
      cat("Training model...\n")
      train_data <- cbind(train_features, response = train_response)
      res <- randomForest(response ~ .,
                          data = train_data,
                          mtry = round((ncol(train_data)-2)/3)) # train each tree
                                                                # on 1/3 features
      if(save_training_data) {
        saveRDS(list(result = res,
                     train_features = train_features,
                     train_response = train_response,
                     test_features = test_features,
                     test_response = test_response), save_fn)
      } else {
        saveRDS(list(result = res), save_fn)
      }
    }

    if(output_weights) {
      if(file.exists(save_fn)) {
        res_obj <- readRDS(save_fn)
        res <- res_obj$result
        train_features <- res_obj$train_features
        train_response <- res_obj$train_response
        test_features <- res_obj$test_features
        test_response <- res_obj$test_response
      } else {
        stop(paste0("Model not found at: ", save_fn, "!\n"))
      }
      saveRDS(varImp(res),
              file.path(save_dir,
                        paste0(use_baseline,
                               "_",
                               use_result_type,
                               ifelse(length(DE_methods) == 1, paste0("_", DE_methods), ""),
                               "_variable-importance.rds")))
    }
  }
}
