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
  
  # COMP_LOG_MEAN_D <- abs(mean(log(counts_A + pseudocount)) - mean(log(counts_B + pseudocount)))
  
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
#' @param use_totals pull fold change in total estimate from absolute counts
#' @param use_renorm_counts pull features generated from method-renormalized
#' data
#' @param use_cpm use counts per million for relative abundances
#' @return data.frame with predictive model training features
#' @import dplyr
#' @import RSQLite
#' @export
pull_features <- function(DE_methods = c("ALDEx2", "DESeq2", "scran"),
                          use_baseline = "oracle",
                          use_totals = FALSE,
                          use_renorm_counts = FALSE,
                          use_cpm = FALSE) {

  conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
  
  # Even though I'm using a SQLite database, it seems much quicker to pull 
  # individual tables and join them in dplyr (ironically) than to do some in the
  # SQL query.
  
  # Note: I'm not pulling FW_RA_PFC1_D, FW_CLR_MED_D, FW_CLR_SD_D, FW_CLR_PNEG_D
  # because these predictors appear to be strongly correlated.
  
  if(use_totals) {
    datasets <- dbGetQuery(conn, paste0("SELECT UUID, P, FC_ABSOLUTE, PERCENT_DIFF_REALIZ FROM datasets"))
  } else {
    datasets <- dbGetQuery(conn, paste0("SELECT UUID, P FROM datasets"))
  }
  if(use_renorm_counts) {
    characteristics <- dbGetQuery(conn, paste0("SELECT * FROM characteristics WHERE PARTIAL=0"))
  } else {
    if(use_cpm) {
      characteristics <- dbGetQuery(conn, paste0("SELECT * FROM characteristics WHERE PARTIAL=0 AND TYPE='cpm'"))
    } else {
      characteristics <- dbGetQuery(conn, paste0("SELECT * FROM characteristics WHERE PARTIAL=0 AND TYPE='relative_abundances'"))
    }
  }
  results <- dbGetQuery(conn, paste0("SELECT UUID, METHOD, TPR, FPR FROM results WHERE PARTIAL_INFO=0 AND BASELINE_TYPE='", use_baseline, "'"))
    
  # Subset to reasonable combinations (i.e. we don't care about joining features
  # for counts that have been scaled by ALDEx2 and had differential abundance
  # called on them by scran because we won't be training the predictive model on
  # that scenario).
  combo_df <- datasets %>%
    left_join(characteristics, by = "UUID") %>%
    left_join(results, by = "UUID") %>%
    filter(METHOD %in% DE_methods) %>%
    select(-c(PARTIAL, FW_RA_PFC1_D, FW_CLR_MED_D, FW_CLR_SD_D, FW_CLR_PNEG_D)) %>%
    filter(TYPE == "relative_abundances" |
             (TYPE == "scaled_ALDEx2" & METHOD == "ALDEx2") |
             (TYPE == "scaled_DESeq2" & METHOD == "DESeq2") |
             (TYPE == "scaled_scran" & METHOD == "scran"))
  
  dbDisconnect(conn)
  
  # Strip "result-less" entries which occasionally result from simulations
  # which don't have any meaningful differential abundance
  combo_df <- combo_df %>%
    filter(!is.na(TPR) & !is.na(FPR))
  
  return(combo_df)
}

#' Generate simulated differential expression for two conditions
#'
#' @param DE_methods which DE calling method's results to predict; if "all",
#' prediction is over all results together
#' @param use_baseline "self" or "oracle"
#' @param use_totals pull fold change in total estimate from absolute counts
#' @param use_renorm_counts pull features generated from method-renormalized
#' data
#' @param output_weights flag indicating whether or not to plot some visualization
#' of feature weights
#' @param train_percent percent of simulated datasets to train on
#' @param use_cpm use counts per million for relative abundances
#' @return NULL (fitted models are saved in output directory)
#' @import randomForest
#' @import dplyr
#' @import tidyr
#' @export
fit_predictive_model <- function(DE_methods = c("ALDEx2", "DESeq2", "scran"),
                                 use_baseline = "oracle",
                                 use_totals = FALSE,
                                 use_renorm_counts = FALSE,
                                 output_weights = TRUE,
                                 train_percent = 0.8,
                                 use_cpm = FALSE) {
  
  if(!(use_baseline %in% c("self", "oracle"))) {
    stop(paste0("Invalid baseline: ", use_baseline, "!\n"))
  }
  
  save_path <- list("output",
                    "predictive_fits",
                    use_baseline,
                    "regression")
  for(i in 1:length(save_path)) {
    fp <- do.call(file.path, save_path[1:i])
    if(!dir.exists(fp)) {
      dir.create(fp)
    }
  }
  save_dir <- fp
  
  features <- pull_features(DE_methods = DE_methods,
                            use_baseline = use_baseline,
                            use_totals = use_totals,
                            use_renorm_counts = use_renorm_counts,
                            use_cpm = use_cpm)

  if(use_renorm_counts) {
    # Combine features from relative abundances and renormalized counts for the 
    # same data set
    f1 <- features %>%
      filter(TYPE == "relative_abundances")
    f2 <- features %>%
      filter(TYPE != "relative_abundances")
    
    # Drop duplicate columns
    f1 <- f1 %>%
      select(-c(P, TYPE, METHOD, TPR, FPR))
    if(use_totals) {
      f1 <- f1 %>%
        select(-c(FC_ABSOLUTE, PERCENT_DIFF_REALIZ))
    }
    
    ignore_columns <- c("UUID", "P", "FC_ABSOLUTE", "PERCENT_DIFF_REALIZ", "TYPE", "METHOD", "TPR", "FPR")
    rename_flag <- !(colnames(f1) %in% ignore_columns)
    colnames(f1)[rename_flag] <- paste0(colnames(f1)[rename_flag], "_rel")
    rename_flag <- !(colnames(f2) %in% ignore_columns)
    colnames(f2)[rename_flag] <- paste0(colnames(f2)[rename_flag], "_scaled")
    
    features <- f1 %>%
      left_join(f2, by = "UUID")
  }
  
  features <- features %>%
    select(-c(TYPE))
  
  if(length(DE_methods) == 1) {
    features <- features %>%
      select(-c(METHOD))
  }
  
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

    features_result_type <- features_result_type %>%
      select(-c(UUID))
    
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
    
    save_fn <- file.path(save_dir,
                         paste0(use_result_type,
                                ifelse(length(DE_methods) == 1, paste0("_", DE_methods), ""),
                                ".rds"))
    
    cat(paste0("Predictive model saving/saved to: ", save_fn, "\n"))
    # if(!file.exists(save_fn)) {
      cat("Training model...\n")
      train_data <- cbind(train_features, response = train_response)
      res <- randomForest(response ~ .,
                          data = train_data)
      # mtry = round((ncol(train_data)-2)/3)) # train each tree
      #                                       # on 1/3 features
      # Save training data
      saveRDS(list(result = res,
                   train_features = train_features,
                   train_response = train_response,
                   test_features = test_features,
                   test_response = test_response), save_fn)
    # }

    if(output_weights) {
      if(file.exists(save_fn)) {
        res_obj <- readRDS(save_fn)
        res <- res_obj$result
        train_features <- res$train_features
        train_response <- res$train_response
        test_features <- res$test_features
        test_response <- res$test_response
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

#' Scale counts using ALDEx2-like size factors
#'
#' @param counts count table oriented as samples (rows) x features (columns)
#' @param pseudocount small count to augment whole matrix by in order to 
#' eliminate zero counts
#' @return scaled count table oriented as samples (rows) x features (columns)
#' @import ALDEx2
#' @export
scaled_counts_ALDEx2 <- function(counts, pseudocount = 0.5) {
  counts <- t(counts)
  scaled_counts <- apply(counts, 2, function(x) {
    # Compute IQR
    bounds <- quantile(x, probs = c(0.25, 0.75))
    # Pull counts/features within these boundaries
    y <- x[x >= bounds[1] & x <= bounds[2]]
    # Compute geometric mean of these counts/features
    sf <- exp(mean(log(y + pseudocount)))
    # Scale counts
    x/sf
  })
  t(scaled_counts)
}

#' Scale counts using DESeq2-like size factors
#'
#' @param counts count table oriented as samples (rows) x features (columns)
#' @param groups group or treatment labels for columns
#' @param pseudocount small count to augment whole matrix by in order to 
#' eliminate zero counts
#' @return scaled count table oriented as samples (rows) x features (columns)
#' @import DESeq2
#' @export
scaled_counts_DESeq2 <- function(counts, groups, pseudocount = 0.5) {
  counts <- t(counts)
  coldata <- data.frame(treatment = groups)
  coldata$treatment <- factor(coldata$treatment)
  levels(coldata$treatment) <- c("control", "treatment")
  dds <- suppressMessages(DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ treatment))
  dds <- estimateSizeFactors(dds)
  sf <- sizeFactors(dds)
  scaled_counts <- t(counts)/sf
  scaled_counts
}

#' Scale counts using scran-like size factors
#'
#' @param counts count table oriented as samples (rows) x features (columns)
#' @param groups group or treatment labels for columns
#' @param pseudocount small count to augment whole matrix by in order to 
#' eliminate zero counts
#' @return scaled count table oriented as samples (rows) x features (columns)
#' @import scran
#' @export
scaled_counts_scran <- function(counts, groups, pseudocount = 0.5) {
  counts <- t(counts)
  count_table <- counts # features x samples
  n_genes <- nrow(count_table)
  cell_metadata <- data.frame(condition = groups)
  rownames(cell_metadata) <- colnames(count_table)
  rownames(count_table) <- paste0("gene", 1:n_genes)
  colnames(count_table) <- paste0("cell", 1:ncol(count_table))
  sce <- SingleCellExperiment(assays = list(counts = count_table),
                              colData = cell_metadata)
  # Specify clusters
  clusters <- as.numeric(as.factor(groups))
  if(min(clusters) == 0) {
    clusters <- clusters + 1
  }
  sce <- suppressWarnings(computeSumFactors(sce, clusters = clusters))
  sf <- sizeFactors(sce)
  scaled_counts <- t(counts)/sf
  scaled_counts
}
