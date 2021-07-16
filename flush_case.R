source("path_fix.R")

library(RSQLite)

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
dbExecute(conn, "DELETE FROM datasets WHERE P=5000 AND CORRP=1")

# Missing in DB
output_dir <- file.path("output", "datasets")
simulation_fn <- list.files(output_dir, pattern = ".rds")
fs_uuids <- unname(sapply(simulation_fn, function(x) substr(x, 1, 36)))
ds_uuids <- dbGetQuery(conn, "SELECT DISTINCT UUID FROM datasets")$UUID
char_uuids <- dbGetQuery(conn, "SELECT DISTINCT UUID FROM characteristics")$UUID
res_uuids <- dbGetQuery(conn, "SELECT DISTINCT UUID FROM results")$UUID
db_uuids <- unique(c(ds_uuids, char_uuids, res_uuids))

unmatched_set <- setdiff(fs_uuids, db_uuids)
cat("Found",length(unmatched_set),"orphaned files...\n")
if(length(unmatched_set) > 0) {
  for(uuid in unmatched_set) {
    unmatched_fn <- file.path(output_dir, paste0(uuid, ".rds"))
    discard <- file.remove(unmatched_fn)
  }
}

dbDisconnect(conn)

# Need to clean up this script.
# Here's the old "orphan" code from a different script.

# simulation_fn <- list.files(output_dir, pattern = ".rds")
# fs_uuids <- unname(sapply(simulation_fn, function(x) substr(x, 1, 36)))
# 
# # Remove DB entries w/ missing files
# ds_uuids <- dbGetQuery(conn, "SELECT DISTINCT UUID FROM datasets")$UUID
# char_uuids <- dbGetQuery(conn, "SELECT DISTINCT UUID FROM characteristics")$UUID
# res_uuids <- dbGetQuery(conn, "SELECT DISTINCT UUID FROM results")$UUID
# db_uuids <- unique(c(ds_uuids, char_uuids, res_uuids))
# 
# # Missing locally
# unmatched_set <- setdiff(db_uuids, fs_uuids)
# cat("Found",length(unmatched_set),"orphaned DB entries...\n")
# if(length(unmatched_set) > 0) {
#   uuid_list_str <- paste0("('",paste(unmatched_set, collapse = "','"),"')")
#   discard <- dbExecute(conn, paste0("DELETE FROM datasets WHERE uuid IN ",
#                                     uuid_list_str,";"))
#   discard <- dbExecute(conn, paste0("DELETE FROM characteristics WHERE uuid IN ",
#                                     uuid_list_str,";"))
#   discard <- dbExecute(conn, paste0("DELETE FROM results WHERE uuid IN ",
#                                     uuid_list_str,";"))
# }
# 
# # Missing in DB
# unmatched_set <- setdiff(fs_uuids, db_uuids)
# cat("Found",length(unmatched_set),"orphaned files...\n")
# if(length(unmatched_set) > 0) {
#   for(uuid in unmatched_set) {
#     unmatched_fn <- file.path(output_dir, paste0(uuid, ".rds"))
#     discard <- file.remove(unmatched_fn)
#   }
# }
