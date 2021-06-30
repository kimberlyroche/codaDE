source("path_fix.R")

library(RSQLite)
library(tidyverse)
library(codaDE)
library(entropy)
library(driver)
library(moments)

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
output_dir <- file.path("output", "datasets")

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

# Update if 1) dataset + either PARTIAL condition is missing an entry in characteristics table or
#           2) if those entries exist but have missing values
# Check PARTIAL = 0 condition
uuid_status <- dbGetQuery(conn, "SELECT * FROM datasets LEFT JOIN characteristics ON datasets.UUID=characteristics.UUID")
update_vector <- apply(uuid_status %>% select(!c(UUID, P, CORRP, TIMESTAMP, PARTIAL)), 1, function(x) sum(is.na(x)) > 0)
update_uuids <- uuid_status$UUID[update_vector]
# Update all
# update_uuids <- dbGetQuery(conn,
#                            paste0("SELECT DISTINCT UUID FROM datasets;"))$UUID

update_uuids <- c("ddae19f2-3ce9-4a0e-b317-3ab8e4dd5c54",
                  "39f1bba3-1d30-43f7-835c-e06c8fb08134")

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
    counts_A <- counts[1:n,]
    counts_B <- counts[(n+1):(n*2),]
    
    features <- characterize_dataset(counts_A, counts_B)

    # --------------------------------------------------------------------------
    #   Record
    # --------------------------------------------------------------------------
    
    # If exists, update
    results <- dbGetQuery(conn,
                          paste0("SELECT * FROM characteristics WHERE UUID='",
                                 uuid, "' AND PARTIAL = ", partial, ";"))
    
    if(nrow(results) > 0) {
      # UPDATE
      discard <- dbExecute(conn, paste0("UPDATE characteristics SET ",
                                        "TOTALS_C_FC=",features$TOTALS_C_FC,", ",
                                        "TOTALS_C_D=",features$TOTALS_C_D,", ",
                                        "TOTALS_C_MAX_D=",features$TOTALS_C_MAX_D,", ",
                                        "TOTALS_C_MED_D=",features$TOTALS_C_MED_D,", ",
                                        "TOTALS_C_SD_D=",features$TOTALS_C_SD_D,", ",
                                        "CORR_RA_MED=",features$CORR_RA_MED,", ",
                                        "CORR_RA_SD=",features$CORR_RA_SD,", ",
                                        "CORR_RA_SKEW=",features$CORR_RA_SKEW,", ",
                                        "CORR_LOG_MED=",features$CORR_LOG_MED,", ",
                                        "CORR_LOG_SD=",features$CORR_LOG_SD,", ",
                                        "CORR_LOG_SKEW=",features$CORR_LOG_SKEW,", ",
                                        "CORR_CLR_MED=",features$CORR_CLR_MED,", ",
                                        "CORR_CLR_SD=",features$CORR_CLR_SD,", ",
                                        "CORR_CLR_SKEW=",features$CORR_CLR_SKEW,", ",
                                        "COMP_C_P0_A=",features$COMP_C_P0_A,", ",
                                        "COMP_C_P0_B=",features$COMP_C_P0_B,", ",
                                        "COMP_C_P1_A=",features$COMP_C_P1_A,", ",
                                        "COMP_C_P1_B=",features$COMP_C_P1_B,", ",
                                        "COMP_C_P5_A=",features$COMP_C_P5_A,", ",
                                        "COMP_C_P5_B=",features$COMP_C_P5_B,", ",
                                        "COMP_RA_P01_A=",features$COMP_RA_P01_A,", ",
                                        "COMP_RA_P01_B=",features$COMP_RA_P01_B,", ",
                                        "COMP_RA_P1_A=",features$COMP_RA_P1_A,", ",
                                        "COMP_RA_P1_B=",features$COMP_RA_P1_B,", ",
                                        "COMP_RA_P5_A=",features$COMP_RA_P5_A,", ",
                                        "COMP_RA_P5_B=",features$COMP_RA_P5_B,", ",
                                        "COMP_RA_MAX_A=",features$COMP_RA_MAX_A,", ",
                                        "COMP_RA_MED_A=",features$COMP_RA_MED_A,", ",
                                        "COMP_RA_SD_A=",features$COMP_RA_SD_A,", ",
                                        "COMP_RA_SKEW_A=",features$COMP_RA_SKEW_A,", ",
                                        "COMP_RA_MAX_B=",features$COMP_RA_MAX_B,", ",
                                        "COMP_RA_MED_B=",features$COMP_RA_MED_B,", ",
                                        "COMP_RA_SD_B=",features$COMP_RA_SD_B,", ",
                                        "COMP_RA_SKEW_B=",features$COMP_RA_SKEW_B,", ",
                                        "COMP_C_ENT_A=",features$COMP_C_ENT_A,", ",
                                        "COMP_C_ENT_B=",features$COMP_C_ENT_B,", ",
                                        "FW_RA_MAX_D=",features$FW_RA_MAX_D,", ",
                                        "FW_RA_MED_D=",features$FW_RA_MED_D,", ",
                                        "FW_RA_SD_D=",features$FW_RA_SD_D,", ",
                                        "FW_RA_PPOS_D=",features$FW_RA_PPOS_D,", ",
                                        "FW_RA_PNEG_D=",features$FW_RA_PNEG_D,", ",
                                        "FW_RA_PFC05_D=",features$FW_RA_PFC05_D,", ",
                                        "FW_RA_PFC1_D=",features$FW_RA_PFC1_D,", ",
                                        "FW_RA_PFC2_D=",features$FW_RA_PFC2_D,", ",
                                        "FW_LOG_MAX_D=",features$FW_LOG_MAX_D,", ",
                                        "FW_LOG_MED_D=",features$FW_LOG_MED_D,", ",
                                        "FW_LOG_SD_D=",features$FW_LOG_SD_D,", ",
                                        "FW_LOG_PPOS_D=",features$FW_LOG_PPOS_D,", ",
                                        "FW_LOG_PNEG_D=",features$FW_LOG_PNEG_D,", ",
                                        "FW_LOG_PFC05_D=",features$FW_LOG_PFC05_D,", ",
                                        "FW_LOG_PFC1_D=",features$FW_LOG_PFC1_D,", ",
                                        "FW_LOG_PFC2_D=",features$FW_LOG_PFC2_D,", ",
                                        "FW_CLR_MAX_D=",features$FW_CLR_MAX_D,", ",
                                        "FW_CLR_MED_D=",features$FW_CLR_MED_D,", ",
                                        "FW_CLR_SD_D=",features$FW_CLR_SD_D,", ",
                                        "FW_CLR_PPOS_D=",features$FW_CLR_PPOS_D,", ",
                                        "FW_CLR_PNEG_D=",features$FW_CLR_PNEG_D,", ",
                                        "FW_CLR_PFC05_D=",features$FW_CLR_PFC05_D,", ",
                                        "FW_CLR_PFC1_D=",features$FW_CLR_PFC1_D,", ",
                                        "FW_CLR_PFC2_D=",features$FW_CLR_PFC2_D," ",
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
                                  "'", uuid, "', ",partial,", ",
                                  features$TOTALS_C_FC, ", ",
                                  features$TOTALS_C_D, ", ",
                                  features$TOTALS_C_MAX_D, ", ",
                                  features$TOTALS_C_MED_D, ", ",
                                  features$TOTALS_C_SD_D, ", ",
                                  features$CORR_RA_MED, ", ",
                                  features$CORR_RA_SD, ", ",
                                  features$CORR_RA_SKEW, ", ",
                                  features$CORR_LOG_MED, ", ",
                                  features$CORR_LOG_SD, ", ",
                                  features$CORR_LOG_SKEW, ", ",
                                  features$CORR_CLR_MED, ", ",
                                  features$CORR_CLR_SD, ", ",
                                  features$CORR_CLR_SKEW, ", ",
                                  features$COMP_C_P0_A, ", ",
                                  features$COMP_C_P0_B, ", ",
                                  features$COMP_C_P1_A, ", ",
                                  features$COMP_C_P1_B, ", ",
                                  features$COMP_C_P5_A, ", ",
                                  features$COMP_C_P5_B, ", ",
                                  features$COMP_RA_P01_A, ", ",
                                  features$COMP_RA_P01_B, ", ",
                                  features$COMP_RA_P1_A, ", ",
                                  features$COMP_RA_P1_B, ", ",
                                  features$COMP_RA_P5_A, ", ",
                                  features$COMP_RA_P5_B, ", ",
                                  features$COMP_RA_MAX_A, ", ",
                                  features$COMP_RA_MED_A, ", ",
                                  features$COMP_RA_SD_A, ", ",
                                  features$COMP_RA_SKEW_A, ", ",
                                  features$COMP_RA_MAX_B, ", ",
                                  features$COMP_RA_MED_B, ", ",
                                  features$COMP_RA_SD_B, ", ",
                                  features$COMP_RA_SKEW_B, ", ",
                                  features$COMP_C_ENT_A, ", ",
                                  features$COMP_C_ENT_B, ", ",
                                  features$FW_RA_MAX_D, ", ",
                                  features$FW_RA_MED_D, ", ",
                                  features$FW_RA_SD_D, ", ",
                                  features$FW_RA_PPOS_D, ", ",
                                  features$FW_RA_PNEG_D, ", ",
                                  features$FW_RA_PFC05_D, ", ",
                                  features$FW_RA_PFC1_D, ", ",
                                  features$FW_RA_PFC2_D, ", ",
                                  features$FW_LOG_MAX_D, ", ",
                                  features$FW_LOG_MED_D, ", ",
                                  features$FW_LOG_SD_D, ", ",
                                  features$FW_LOG_PPOS_D, ", ",
                                  features$FW_LOG_PNEG_D, ", ",
                                  features$FW_LOG_PFC05_D, ", ",
                                  features$FW_LOG_PFC1_D, ", ",
                                  features$FW_LOG_PFC2_D, ", ",
                                  features$FW_CLR_MAX_D, ", ",
                                  features$FW_CLR_MED_D, ", ",
                                  features$FW_CLR_SD_D, ", ",
                                  features$FW_CLR_PPOS_D, ", ",
                                  features$FW_CLR_PNEG_D, ", ",
                                  features$FW_CLR_PFC05_D, ", ",
                                  features$FW_CLR_PFC1_D, ", ",
                                  features$FW_CLR_PFC2_D,
                                  ");"))
    }
  }
}
# dbGetQuery(conn, "SELECT * FROM characteristics;")

dbDisconnect(conn)
