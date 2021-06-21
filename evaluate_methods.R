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

n <- opt$n
p <- opt$p
method <- opt$method

if(!(method %in% allowed_methods)) {
  stop(paste0("Method '",method,"' not allowed!\n"))
}

# Randomly wait a few seconds before connecting. This is a stupid but effective
# way to prevent a race condition that seems to allow the first job in a batch
# to connect and write, but not subsequent jobs in a batch.

Sys.sleep(runif(1, min = 1, max = 10))

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
# Increase the "busy" timeout; default is too short
discard <- dbExecute(conn, "PRAGMA busy_timeout = 60000;")

output_dir <- file.path("output", "datasets")

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

cat(paste0("Found ", length(all_uuids), " datasets matching criteria!\n"))

datasets_to_process <- data.frame(uuid = c(),
                                  run_id = c(),
                                  result_type = c(),
                                  partial_info = c())
for(uuid in sample(all_uuids)) {
  run_id <- dbGetQuery(conn, paste0("SELECT RUN_ID FROM results ",
                                    "WHERE UUID = '",uuid,"' AND ",
                                    "METHOD = '",method,"' AND ",
                                    "(RESULT IS NOT NULL OR FIT_IN_PROGRESS = 1);"))
  if(nrow(run_id) == 0) {
    cat(paste0("Selecting candidate ", uuid, "\n"))
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
  } else {
    cat(paste0("Candidate ", uuid, " has already been processed!\n"))
  }
  if(nrow(datasets_to_process) > 0 && datasets_to_process %>% distinct(uuid) %>% tally() %>% pull(n) == n) {
    break;
  }
}
discard <- dbExecute(conn, "COMMIT;")

for(this_uuid in unique(datasets_to_process$uuid)) {
  cat("Updating UUID: ", this_uuid, "\n")
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
  discard <- dbExecute(conn, paste0("UPDATE results SET RESULT = ",rates$tpr,", ",
                                    "FIT_IN_PROGRESS = 0 ",
                                    "WHERE RUN_ID = ",tpr_run_id,";"))
  cat(paste0("Updated ", discard, " entry with run ID ", tpr_run_id, "\n"))
  discard <- dbExecute(conn, paste0("UPDATE results SET RESULT = ",rates$fpr,", ",
                                    "FIT_IN_PROGRESS = 0 ",
                                    "WHERE RUN_ID = ",fpr_run_id,";"))
  cat(paste0("Updated ", discard, " entry with run ID ", fpr_run_id, "\n"))
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
  discard <- dbExecute(conn, paste0("UPDATE results SET RESULT = ",rates$tpr,", ",
                                    "FIT_IN_PROGRESS = 0 ",
                                    "WHERE RUN_ID = ",tpr_partial_run_id,";"))
  cat(paste0("Updated ", discard, " entry with run ID ", tpr_partial_run_id, "\n"))
  discard <- dbExecute(conn, paste0("UPDATE results SET RESULT = ",rates$fpr,", ",
                                    "FIT_IN_PROGRESS = 0 ",
                                    "WHERE RUN_ID = ",fpr_partial_run_id,";"))
  cat(paste0("Updated ", discard, " entry with run ID ", fpr_partial_run_id, "\n"))
  discard <- dbExecute(conn, "COMMIT;")
}

# dbGetQuery(conn, "SELECT * FROM results;")
# print(dbGetQuery(conn, paste0("SELECT * FROM results WHERE RUN_ID >= ", tpr_run_id, " AND RUN_ID <= ", fpr_partial_run_id, ";")))

dbDisconnect(conn)
