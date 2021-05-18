source("path_fix.R")

library(RSQLite)
library(tidyverse)
library(codaDE)

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
output_dir <- file.path("output", "datasets")

# ------------------------------------------------------------------------------
#   Sync folder and DB
# ------------------------------------------------------------------------------

simulation_fn <- list.files(output_dir, pattern = ".rds")
fs_uuids <- unname(sapply(simulation_fn, function(x) substr(x, 1, 36)))

# Remove DB entries w/ missing files
db_uuids <- dbGetQuery(conn, "SELECT UUID FROM datasets")$UUID

# Missing locally
unmatched_set <- setdiff(db_uuids, fs_uuids)
cat("Found",length(unmatched_set),"orphaned DB entries...\n")
if(length(unmatched_set) > 0) {
  uuid_list_str <- paste0("('",paste(unmatched_set, collapse = "','"),"')")
  # dbGetQuery(conn, paste0("SELECT * FROM datasets WHERE uuid IN ",uuid_list_str,";"))
  discard <- dbExecute(conn, paste0("DELETE FROM datasets WHERE uuid IN ",uuid_list_str,";"))
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

# Characterize data: add observed fold change, mean and median observed correlation
update_uuids <- dbGetQuery(conn, paste0("SELECT UUID FROM datasets WHERE ",
                                        "FOLD_CHANGE IS NULL OR ",
                                        "FOLD_CHANGE_PARTIAL IS NULL OR ",
                                        "MEAN_CORR IS NULL OR ",
                                        "MEAN_CORR_PARTIAL IS NULL OR ",
                                        "MEDIAN_CORR IS NULL OR ",
                                        "MEDIAN_CORR_PARTIAL IS NULL;"))$UUID
for(uuid in update_uuids) {
  fn <- file.path(output_dir, paste0(uuid, ".rds"))
  dataset <- readRDS(fn)
  # Calculate fold change - no/partial information cases
  counts1 <- dataset$simulation$observed_counts1
  counts2 <- dataset$simulation$observed_counts2
  half_n <- nrow(counts1)/2
  m1 <- mean(rowSums(counts1[1:n,]))
  m2 <- mean(rowSums(counts1[(n+1):(n*2),]))
  fold_diff1 <- max(c(m1, m2)) / min(c(m1, m2))
  m1 <- mean(rowSums(counts2[1:n,]))
  m2 <- mean(rowSums(counts2[(n+1):(n*2),]))
  fold_diff2 <- max(c(m1, m2)) / min(c(m1, m2))
  # Calculate observed correlation between features - no/partial information cases
  observed_correlation <- cor(log(counts1 + 0.5))
  observed_correlation <- observed_correlation[upper.tri(observed_correlation,
                                                         diag = FALSE)]
  mean_cor1 <- mean(observed_correlation)
  median_cor1 <- median(observed_correlation)
  observed_correlation <- cor(log(counts2 + 0.5))
  observed_correlation <- observed_correlation[upper.tri(observed_correlation,
                                                         diag = FALSE)]
  mean_cor2 <- mean(observed_correlation)
  median_cor2 <- median(observed_correlation)
  
  # Update DB
  discard <- dbExecute(conn, paste0("UPDATE datasets SET ",
                                    "FOLD_CHANGE = ",fold_diff1,", ",
                                    "FOLD_CHANGE_PARTIAL = ",fold_diff2,", ",
                                    "MEAN_CORR = ",mean_cor1,", ",
                                    "MEAN_CORR_PARTIAL = ",mean_cor2,", ",
                                    "MEDIAN_CORR = ",median_cor1,", ",
                                    "MEDIAN_CORR_PARTIAL = ",median_cor2," ",
                                    "WHERE UUID == '", uuid, "';"))
}

dbGetQuery(conn, "SELECT * FROM datasets;")

dbDisconnect(conn)

