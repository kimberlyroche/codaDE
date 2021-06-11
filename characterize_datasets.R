source("path_fix.R")

library(RSQLite)
library(tidyverse)
library(codaDE)
library(entropy)
library(driver)
library(moments)

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
output_dir <- file.path("output", "datasets")
pseudocount <- 0.1

# ------------------------------------------------------------------------------
#   Sync folder and DB
# ------------------------------------------------------------------------------

simulation_fn <- list.files(output_dir, pattern = ".rds")
fs_uuids <- unname(sapply(simulation_fn, function(x) substr(x, 1, 36)))

# Remove DB entries w/ missing files
ds_uuids <- dbGetQuery(conn, "SELECT DISTINCT UUID FROM datasets")$UUID
char_uuids <- dbGetQuery(conn, "SELECT DISTINCT UUID FROM characteristics")$UUID
res_uuids <- dbGetQuery(conn, "SELECT DISTINCT UUID FROM results")$UUID
db_uuids <- unique(c(ds_uuids, char_uuids, res_uuids))

# Missing locally
unmatched_set <- setdiff(db_uuids, fs_uuids)
cat("Found",length(unmatched_set),"orphaned DB entries...\n")
if(length(unmatched_set) > 0) {
  uuid_list_str <- paste0("('",paste(unmatched_set, collapse = "','"),"')")
  discard <- dbExecute(conn, paste0("DELETE FROM datasets WHERE uuid IN ",
                                    uuid_list_str,";"))
  discard <- dbExecute(conn, paste0("DELETE FROM characteristics WHERE uuid IN ",
                                    uuid_list_str,";"))
  discard <- dbExecute(conn, paste0("DELETE FROM results WHERE uuid IN ",
                                    uuid_list_str,";"))
}

# Missing in DB
unmatched_set <- setdiff(fs_uuids, db_uuids)
cat("Found",length(unmatched_set),"orphaned files...\n")
if(length(unmatched_set) > 0) {
  for(uuid in unmatched_set) {
    unmatched_fn <- file.path(output_dir, paste0(uuid, ".rds"))
    discard <- file.remove(unmatched_fn)
  }
}

# If any features are missing, we'll update them all
update_uuids <- dbGetQuery(conn,
                           paste0("SELECT DISTINCT UUID FROM characteristics ",
                                  "WHERE ",
                                  "FOLD_CHANGE IS NULL OR ",
                                  "MEAN_CORR IS NULL OR ",
                                  "MEDIAN_CORR IS NULL OR ",
                                  "BASE_SPARSITY IS NULL OR ",
                                  "DELTA_SPARSITY IS NULL OR ",
                                  "PERCENT_STABLE IS NULL OR ",
                                  "DIR_CONSENSUS IS NULL OR ",
                                  "MAX_DELTA_REL IS NULL OR ",
                                  "MEDIAN_DELTA_REL IS NULL OR ",
                                  "BASE_ENTROPY IS NULL OR ",
                                  "DELTA_ENTROPY IS NULL;"))$UUID
for(u in 1:length(update_uuids)) {
  uuid <- update_uuids[u]
  cat(paste0("Updating ", uuid, " (", u, "/", length(update_uuids), ")\n"))
  fn <- file.path(output_dir, paste0(uuid, ".rds"))
  dataset <- readRDS(fn)
  
  for(partial in c(0, 1)) {
    if(partial == 1) {
      counts <- dataset$simulation$observed_counts2
    } else {
      counts <- dataset$simulation$observed_counts1
    }
    n <- nrow(counts)/2
    p <- ncol(counts)
    counts_A <- counts[1:n,]
    counts_B <- counts[(n+1):(n*2),]
    
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
    
    COMP_C_P0_A <- sum(counts_A == 0) / (n*p)
    COMP_C_P0_B <- sum(counts_B == 0) / (n*p)
    COMP_C_P1_A <- sum(counts_A == 1) / (n*p)
    COMP_C_P1_B <- sum(counts_B == 1) / (n*p)
    COMP_C_P5_A <- sum(counts_A <= 5) / (n*p)
    COMP_C_P5_B <- sum(counts_B <= 5) / (n*p)
    
    COMP_RA_P01_A <- sum(relab_A < 0.001) / (n*p)
    COMP_RA_P01_B <- sum(relab_B < 0.001) / (n*p)
    COMP_RA_P1_A <- sum(relab_A < 0.01) / (n*p)
    COMP_RA_P1_B <- sum(relab_B < 0.01) / (n*p)
    COMP_RA_P5_A <- sum(relab_A < 0.05) / (n*p)
    COMP_RA_P5_B <- sum(relab_B < 0.05) / (n*p)
    
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
    
    # --------------------------------------------------------------------------
    #   Record
    # --------------------------------------------------------------------------
    
    # If exists, update
    results <- dbGetQuery(conn,
                          paste0("SELECT * FROM characteristics WHERE UUID='",
                                 uuid, "' AND PARTIAL = 0;"))
    
    if(nrow(results) > 0) {
      # UPDATE
      discard <- dbExecute(conn, paste0("UPDATE characteristics SET ",
                                        "TOTALS_C_FC=",TOTALS_C_FC,", ",
                                        "TOTALS_C_D=",TOTALS_C_D,", ",
                                        "TOTALS_C_MAX_D=",TOTALS_C_MAX_D,", ",
                                        "TOTALS_C_MED_D=",TOTALS_C_MED_D,", ",
                                        "TOTALS_C_SD_D=",TOTALS_C_SD_D,", ",
                                        "CORR_RA_MED=",CORR_RA_MED,", ",
                                        "CORR_RA_SD=",CORR_RA_SD,", ",
                                        "CORR_RA_SKEW=",CORR_RA_SKEW,", ",
                                        "CORR_LOG_MED=",CORR_LOG_MED,", ",
                                        "CORR_LOG_SD=",CORR_LOG_SD,", ",
                                        "CORR_LOG_SKEW=",CORR_LOG_SKEW,", ",
                                        "CORR_CLR_MED=",CORR_CLR_MED,", ",
                                        "CORR_CLR_SD=",CORR_CLR_SD,", ",
                                        "CORR_CLR_SKEW=",CORR_CLR_SKEW,", ",
                                        "COMP_C_P0_A=",COMP_C_P0_A,", ",
                                        "COMP_C_P0_B=",COMP_C_P0_B,", ",
                                        "COMP_C_P1_A=",COMP_C_P1_A,", ",
                                        "COMP_C_P1_B=",COMP_C_P1_B,", ",
                                        "COMP_C_P5_A=",COMP_C_P5_A,", ",
                                        "COMP_C_P5_B=",COMP_C_P5_B,", ",
                                        "COMP_RA_P01_A=",COMP_RA_P01_A,", ",
                                        "COMP_RA_P01_B=",COMP_RA_P01_B,", ",
                                        "COMP_RA_P1_A=",COMP_RA_P1_A,", ",
                                        "COMP_RA_P1_B=",COMP_RA_P1_B,", ",
                                        "COMP_RA_P5_A=",COMP_RA_P5_A,", ",
                                        "COMP_RA_P5_B=",COMP_RA_P5_B,", ",
                                        "COMP_RA_MAX_A=",COMP_RA_MAX_A,", ",
                                        "COMP_RA_MED_A=",COMP_RA_MED_A,", ",
                                        "COMP_RA_SD_A=",COMP_RA_SD_A,", ",
                                        "COMP_RA_SKEW_A=",COMP_RA_SKEW_A,", ",
                                        "COMP_RA_MAX_B=",COMP_RA_MAX_B,", ",
                                        "COMP_RA_MED_B=",COMP_RA_MED_B,", ",
                                        "COMP_RA_SD_B=",COMP_RA_SD_B,", ",
                                        "COMP_RA_SKEW_B=",COMP_RA_SKEW_B,", ",
                                        "COMP_C_ENT_A=",COMP_C_ENT_A,", ",
                                        "COMP_C_ENT_B=",COMP_C_ENT_B,", ",
                                        "FW_RA_MAX_D=",FW_RA_MAX_D,", ",
                                        "FW_RA_MED_D=",FW_RA_MED_D,", ",
                                        "FW_RA_SD_D=",FW_RA_SD_D,", ",
                                        "FW_RA_PPOS_D=",FW_RA_PPOS_D,", ",
                                        "FW_RA_PNEG_D=",FW_RA_PNEG_D,", ",
                                        "FW_RA_PFC05=",FW_RA_PFC05,", ",
                                        "FW_RA_PFC1=",FW_RA_PFC1,", ",
                                        "FW_RA_PFC2=",FW_RA_PFC2,", ",
                                        "FW_LOG_MAX_D=",FW_LOG_MAX_D,", ",
                                        "FW_LOG_MED_D=",FW_LOG_MED_D,", ",
                                        "FW_LOG_SD_D=",FW_LOG_SD_D,", ",
                                        "FW_LOG_PPOS_D=",FW_LOG_PPOS_D,", ",
                                        "FW_LOG_PNEG_D=",FW_LOG_PNEG_D,", ",
                                        "FW_LOG_PFC05_D=",FW_LOG_PFC05_D,", ",
                                        "FW_LOG_PFC1_D=",FW_LOG_PFC1_D,", ",
                                        "FW_LOG_PFC2_D=",FW_LOG_PFC2_D,", ",
                                        "FW_CLR_MAX_D=",FW_CLR_MAX_D,", ",
                                        "FW_CLR_MED_D=",FW_CLR_MED_D,", ",
                                        "FW_CLR_SD_D=",FW_CLR_SD_D,", ",
                                        "FW_CLR_PPOS_D=",FW_CLR_PPOS_D,", ",
                                        "FW_CLR_PNEG_D=",FW_CLR_PNEG_D,", ",
                                        "FW_CLR_PFC05_D=",FW_CLR_PFC05_D,", ",
                                        "FW_CLR_PFC1_D=",FW_CLR_PFC1_D,", ",
                                        "FW_CLR_PFC2_D=",FW_CLR_PFC2_D,", ",
                                        " WHERE ",
                                        "UUID='",uuid,"' AND ",
                                        "PARTIAL=",partial,";"))
    } else {
      # INSERT
      discard <- dbExecute(conn,
                           paste0("INSERT INTO characteristics(UUID, PARTIAL, ",
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
                                  "FW_CLR_PFC2_D) ",
                                  "VALUES(",
                                  "'", uuid, "', 0, ",
                                  TOTALS_C_FC, ", ",
                                  TOTALS_C_D, ", ",
                                  TOTALS_C_MAX_D, ", ",
                                  TOTALS_C_MED_D, ", ",
                                  TOTALS_C_SD_D, ", ",
                                  CORR_RA_MED, ", ",
                                  CORR_RA_SD, ", ",
                                  CORR_RA_SKEW, ", ",
                                  CORR_LOG_MED, ", ",
                                  CORR_LOG_SD, ", ",
                                  CORR_LOG_SKEW, ", ",
                                  CORR_CLR_MED, ", ",
                                  CORR_CLR_SD, ", ",
                                  CORR_CLR_SKEW, ", ",
                                  COMP_C_P0_A, ", ",
                                  COMP_C_P0_B, ", ",
                                  COMP_C_P1_A, ", ",
                                  COMP_C_P1_B, ", ",
                                  COMP_C_P5_A, ", ",
                                  COMP_C_P5_B, ", ",
                                  COMP_RA_P01_A, ", ",
                                  COMP_RA_P01_B, ", ",
                                  COMP_RA_P1_A, ", ",
                                  COMP_RA_P1_B, ", ",
                                  COMP_RA_P5_A, ", ",
                                  COMP_RA_P5_B, ", ",
                                  COMP_RA_MAX_A, ", ",
                                  COMP_RA_MED_A, ", ",
                                  COMP_RA_SD_A, ", ",
                                  COMP_RA_SKEW_A, ", ",
                                  COMP_RA_MAX_B, ", ",
                                  COMP_RA_MED_B, ", ",
                                  COMP_RA_SD_B, ", ",
                                  COMP_RA_SKEW_B, ", ",
                                  COMP_C_ENT_A, ", ",
                                  COMP_C_ENT_B, ", ",
                                  FW_RA_MAX_D, ", ",
                                  FW_RA_MED_D, ", ",
                                  FW_RA_SD_D, ", ",
                                  FW_RA_PPOS_D, ", ",
                                  FW_RA_PNEG_D, ", ",
                                  FW_RA_PFC05_D, ", ",
                                  FW_RA_PFC1_D, ", ",
                                  FW_RA_PFC2_D, ", ",
                                  FW_LOG_MAX_D, ", ",
                                  FW_LOG_MED_D, ", ",
                                  FW_LOG_SD_D, ", ",
                                  FW_LOG_PPOS_D, ", ",
                                  FW_LOG_PNEG_D, ", ",
                                  FW_LOG_PFC05_D, ", ",
                                  FW_LOG_PFC1_D, ", ",
                                  FW_LOG_PFC2_D, ", ",
                                  FW_CLR_MAX_D, ", ",
                                  FW_CLR_MED_D, ", ",
                                  FW_CLR_SD_D, ", ",
                                  FW_CLR_PPOS_D, ", ",
                                  FW_CLR_PNEG_D, ", ",
                                  FW_CLR_PFC05_D, ", ",
                                  FW_CLR_PFC1_D, ", ",
                                  FW_CLR_PFC2_D,
                                  ");"))
    }
  }
}
# dbGetQuery(conn, "SELECT * FROM characteristics;")

dbDisconnect(conn)
