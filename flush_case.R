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
