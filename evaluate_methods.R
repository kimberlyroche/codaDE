source("path_fix.R")

library(RSQLite)
library(tidyverse)
library(codaDE)
library(optparse)

option_list = list(
  make_option(c("--n"), type = "numeric", default = 10,
              help = "max number of data sets to evaluate", metavar = "numeric"),
  make_option(c("--p"), type = "numeric", default = 0,
              help = "number of genes", metavar = "numeric"),
  make_option(c("--method"), type = "character", default = NULL,
              help = "differential abundance method", metavar = "numeric")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

allowed_methods <- c("NBGLM", "DESeq2", "MAST", "ALDEx2", "scran")

# ------------------------------------------------------------------------------
#  Global params
# ------------------------------------------------------------------------------

# p <- opt$p
# n <- opt$n
# method <- opt$method
p <- 100
n <- 2
method <- "MAST"

if(!(method %in% allowed_methods)) {
  stop(paste0("Method '",method,"' not allowed!\n"))
}

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
output_dir <- file.path("output", "datasets")

# ------------------------------------------------------------------------------
#   Functions
# ------------------------------------------------------------------------------

DE_by_NB <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_NB(ref_data, groups)
  }
  calls <- call_DA_NB(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

DE_by_DESeq2 <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_DESeq2(ref_data, groups)
  }
  calls <- call_DA_DESeq2(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

DE_by_MAST <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_MAST(ref_data, groups)
  }
  calls <- call_DA_MAST(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

DE_by_ALDEx2 <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_ALDEx2(ref_data, groups)
  }
  calls <- call_DA_ALDEx2(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

DE_by_scran <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_scran(ref_data, groups)
  }
  calls <- call_DA_scran(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

calc_DE_discrepancy <- function(ref_data, data, groups, method = "NBGLM",
                                oracle_calls = NULL) {
  if(method == "NBGLM") {
    DE_calls <- DE_by_NB(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
  if(method == "DESeq2") {
    DE_calls <- DE_by_DESeq2(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
  if(method == "MAST") {
    DE_calls <- DE_by_MAST(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
  if(method == "ALDEx2") {
    DE_calls <- DE_by_ALDEx2(ref_data, data, groups,
                             oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
  if(method == "scran") {
    DE_calls <- DE_by_scran(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
  
  de <- calls < 0.05
  sim_de <- oracle_calls < 0.05
  
  TP <- sum(de & sim_de)
  FP <- sum(de & !sim_de)
  
  TN <- sum(!de & !sim_de)
  FN <- sum(!de & sim_de)
  
  return(list(tpr = TP/(TP+FN), fpr = FP/(FP+TN), oracle_calls = oracle_calls))
}

# ------------------------------------------------------------------------------
#   Main logic
# ------------------------------------------------------------------------------

discard <- dbExecute(conn, "BEGIN TRANSACTION;")
if(p < 1) {
  all_uuids <- dbGetQuery(conn, "SELECT UUID FROM datasets;")$UUID
} else {
  all_uuids <- dbGetQuery(conn, paste0("SELECT UUID FROM datasets ",
                                       "WHERE P = ",p,";"))$UUID
}
datasets_to_process <- data.frame(uuid = c(),
                                  run_id = c(),
                                  result_type = c(),
                                  partial_info = c())
for(uuid in sample(all_uuids)) {
  run_id <- dbGetQuery(conn, paste0("SELECT RUN_ID FROM results ",
                                    "WHERE UUID = '",uuid,"' AND ",
                                    "METHOD = '",method,"';"))
  if(nrow(run_id) == 0) {
    # No partial total abundance info case
    discard <- dbExecute(conn,
                         paste0("INSERT INTO results(UUID, METHOD, ",
                                "PARTIAL_INFO, RESULT_TYPE, FIT_IN_PROGRESS) ",
                                "VALUES('",uuid,"','",method,"',0,'tpr',1);"))
    run_id <- dbGetQuery(conn, "SELECT last_insert_rowid() AS RUN_ID;")$RUN_ID
    datasets_to_process <- bind_rows(datasets_to_process,
                                     data.frame(uuid = uuid,
                                                run_id = run_id,
                                                result_type = 'tpr',
                                                partial_info = 0))
    discard <- dbExecute(conn,
                         paste0("INSERT INTO results(UUID, METHOD, ",
                                "PARTIAL_INFO, RESULT_TYPE, FIT_IN_PROGRESS) ",
                                "VALUES('",uuid,"','",method,"',0,'fpr',1);"))
    run_id <- dbGetQuery(conn, "SELECT last_insert_rowid() AS RUN_ID;")$RUN_ID
    datasets_to_process <- bind_rows(datasets_to_process,
                                     data.frame(uuid = uuid,
                                                run_id = run_id,
                                                result_type = 'fpr',
                                                partial_info = 0))
    # Partial total abundance info case
    discard <- dbExecute(conn,
                         paste0("INSERT INTO results(UUID, METHOD, ",
                                "PARTIAL_INFO, RESULT_TYPE, FIT_IN_PROGRESS) ",
                                "VALUES('",uuid,"','",method,"',1,'tpr',1);"))
    run_id <- dbGetQuery(conn, "SELECT last_insert_rowid() AS RUN_ID;")$RUN_ID
    datasets_to_process <- bind_rows(datasets_to_process,
                                     data.frame(uuid = uuid,
                                                run_id = run_id,
                                                result_type = 'tpr',
                                                partial_info = 1))
    discard <- dbExecute(conn,
                         paste0("INSERT INTO results(UUID, METHOD, ",
                                "PARTIAL_INFO, RESULT_TYPE, FIT_IN_PROGRESS) ",
                                "VALUES('",uuid,"','",method,"',1,'fpr',1);"))
    run_id <- dbGetQuery(conn, "SELECT last_insert_rowid() AS RUN_ID;")$RUN_ID
    datasets_to_process <- bind_rows(datasets_to_process,
                                     data.frame(uuid = uuid,
                                                run_id = run_id,
                                                result_type = 'fpr',
                                                partial_info = 1))
  }
  if(nrow(datasets_to_process) > 0 && datasets_to_process %>% distinct(uuid) %>% tally() %>% pull(n) == n) {
    break;
  }
}
discard <- dbExecute(conn, "COMMIT;")

for(this_uuid in unique(datasets_to_process$uuid)) {
  # Get run ID assignments ahead of time
  tpr_run_id <- datasets_to_process %>%
    filter(result_type == 'tpr' & partial_info == 0 & uuid == this_uuid) %>%
    pull(run_id)
  tpr_partial_run_id <- datasets_to_process %>%
    filter(result_type == 'tpr' & partial_info == 1 & uuid == this_uuid) %>%
    pull(run_id)
  fpr_run_id <- datasets_to_process %>%
    filter(result_type == 'fpr' & partial_info == 0 & uuid == this_uuid) %>%
    pull(run_id)
  fpr_partial_run_id <- datasets_to_process %>%
    filter(result_type == 'fpr' & partial_info == 1 & uuid == this_uuid) %>%
    pull(run_id)

  # Parse data set
  data <- readRDS(file.path(output_dir, paste0(this_uuid, ".rds")))
  
  # No partial total abundance information case
  # Call DA using baseline model (NB GLM)
  rates_baseline <- calc_DE_discrepancy(data$simulation$abundances[,1:p],
                                        data$simulation$observed_counts1[,1:p],
                                        data$simulation$groups)
  oracle_calls <- rates_baseline$oracle_calls
  # Call DA using chosen model
  rates <- calc_DE_discrepancy(data$simulation$abundances[,1:p],
                               data$simulation$observed_counts1[,1:p],
                               data$simulation$groups,
                               method = method,
                               oracle_calls = oracle_calls)
  
  discard <- dbExecute(conn, "BEGIN TRANSACTION;")
  discard <- dbExecute(conn, paste0("UPDATE results SET RESULT = ",rates$tpr," ",
                                    "WHERE RUN_ID = ",tpr_run_id,";"))
  discard <- dbExecute(conn, paste0("UPDATE results SET RESULT = ",rates$fpr," ",
                                    "WHERE RUN_ID = ",fpr_run_id,";"))
  discard <- dbExecute(conn, "COMMIT;")
  
  # Partial total abundance information case
  # Call DA using baseline model (NB GLM)
  rates_baseline <- calc_DE_discrepancy(data$simulation$abundances[,1:p],
                                        data$simulation$observed_counts2[,1:p],
                                        data$simulation$groups)
  oracle_calls <- rates_baseline$oracle_calls
  # Call DA using chosen model
  rates <- calc_DE_discrepancy(data$simulation$abundances[,1:p],
                               data$simulation$observed_counts2[,1:p],
                               data$simulation$groups,
                               method = method,
                               oracle_calls = oracle_calls)
  
  discard <- dbExecute(conn, "BEGIN TRANSACTION;")
  discard <- dbExecute(conn, paste0("UPDATE results SET RESULT = ",rates$tpr," ",
                                    "WHERE RUN_ID = ",tpr_partial_run_id,";"))
  discard <- dbExecute(conn, paste0("UPDATE results SET RESULT = ",rates$fpr," ",
                                    "WHERE RUN_ID = ",fpr_partial_run_id,";"))
  discard <- dbExecute(conn, "COMMIT;")
}
# dbGetQuery(conn, "SELECT * FROM results;")

dbDisconnect(conn)
