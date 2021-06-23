source("path_fix.R")

library(codaDE)
library(RSQLite)

# ------------------------------------------------------------------------------
#   Open new and old connections
# ------------------------------------------------------------------------------

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
conn_new <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations_new.db"))

# ------------------------------------------------------------------------------
#   Create tables in new DB
# ------------------------------------------------------------------------------

discard <- dbExecute(conn_new, paste0("CREATE TABLE datasets(",
                                  "UUID VARCHAR(36) PRIMARY KEY,",
                                  "P INTEGER,",
                                  "CORRP INTEGER,",
                                  "TIMESTAMP VARCHAR(19));"))

discard <- dbExecute(conn_new, paste0("CREATE TABLE characteristics(",
                                  "UUID VARCHAR(36),",
                                  "PARTIAL INTEGER,",
                                  "FOLD_CHANGE REAL,",
                                  "MEAN_CORR REAL,",
                                  "MEDIAN_CORR REAL,",
                                  "BASE_SPARSITY REAL,",
                                  "DELTA_SPARSITY REAL,",
                                  "PERCENT_STABLE REAL,",
                                  "DIR_CONSENSUS REAL,",
                                  "MAX_DELTA_REL REAL,",
                                  "MEDIAN_DELTA_REL REAL,",
                                  "BASE_ENTROPY REAL,",
                                  "DELTA_ENTROPY REAL,",
                                  "PRIMARY KEY(UUID, PARTIAL));"))

discard <- dbExecute(conn_new, paste0("CREATE TABLE results(",
                                  "RUN_ID INTEGER PRIMARY KEY AUTOINCREMENT,",
                                  "UUID VARCHAR(36),",
                                  "METHOD VARCHAR(64),",
                                  "PARTIAL_INFO INT,",
                                  "BASELINE VARCHAR(64),",
                                  "RESULT REAL,",
                                  "RESULT_TYPE VARCHAR(24),",
                                  "FIT_IN_PROGRESS INTEGER);"))

# ------------------------------------------------------------------------------
#   Migrate content
# ------------------------------------------------------------------------------

cat("Copying results from OLD datasets table to NEW datasets table...\n")
results <- dbGetQuery(conn, "SELECT UUID, P, CORRP, TIMESTAMP FROM datasets;")
for(i in 1:nrow(results)) {
  line <- results[i,]
  discard <- dbExecute(conn_new,
                       paste0("INSERT INTO datasets(UUID, P, CORRP, TIMESTAMP)",
                              " VALUES('",
                              line$UUID,"', ",
                              line$P,", ",
                              line$CORRP,", '",
                              line$TIMESTAMP,"');"))
}
results <- dbGetQuery(conn, "SELECT COUNT(DISTINCT(UUID)) FROM datasets;")
cat(paste0("Updated ", results[1,1], " UUIDS in datasets table\n"))

cat("Copying results from OLD results table to NEW results table...\n")
# Make sure to omit partially completed results here
results <- dbGetQuery(conn, "SELECT * FROM results WHERE RESULT IS NOT NULL;")
for(i in 1:nrow(results)) {
  line <- results[i,]
  discard <- dbExecute(conn_new,
                       paste0("INSERT INTO results(RUN_ID, UUID, METHOD, ",
                              "PARTIAL_INFO, RESULT, RESULT_TYPE, ",
                              "FIT_IN_PROGRESS) ",
                              "VALUES(",
                              line$RUN_ID,", ",
                              "'",line$UUID,"', ",
                              "'",line$METHOD,"', ",
                              line$PARTIAL_INFO,", ",
                              line$RESULT,", ",
                              "'",line$RESULT_TYPE,"', 0);"))
}
results <- dbGetQuery(conn, "SELECT COUNT(DISTINCT(UUID)) FROM results;")
cat(paste0("Updated ", results[1,1], " UUIDS in results table\n"))

cat("Copying results from OLD datasets table to NEW characteristics table...\n")
# Make sure to omit partially completed results here
results <- dbGetQuery(conn, paste0("SELECT UUID, FOLD_CHANGE, ",
                                   "FOLD_CHANGE_PARTIAL, MEAN_CORR, ",
                                   "MEAN_CORR_PARTIAL, MEDIAN_CORR, ",
                                   "MEDIAN_CORR_PARTIAL FROM datasets WHERE ",
                                   "FOLD_CHANGE IS NOT NULL AND ",
                                   "FOLD_CHANGE_PARTIAL IS NOT NULL AND ",
                                   "MEAN_CORR IS NOT NULL AND ",
                                   "MEAN_CORR_PARTIAL IS NOT NULL AND ",
                                   "MEDIAN_CORR IS NOT NULL AND ",
                                   "MEDIAN_CORR_PARTIAL IS NOT NULL;"))
results_partial <- dbGetQuery(conn, paste0("SELECT UUID FROM datasets WHERE ",
                                   "FOLD_CHANGE IS NULL OR ",
                                   "FOLD_CHANGE_PARTIAL IS NULL OR ",
                                   "FOLD_CHANGE IS NULL OR ",
                                   "MEAN_CORR_PARTIAL IS NULL OR ",
                                   "MEDIAN_CORR IS NULL OR ",
                                   "MEDIAN_CORR_PARTIAL IS NULL;"))
for(i in 1:nrow(results)) {
  line <- results[i,]
  discard <- dbExecute(conn_new, paste0("INSERT INTO characteristics(",
                                        "UUID, PARTIAL, FOLD_CHANGE, ",
                                        "MEAN_CORR, MEDIAN_CORR) ",
                                        "VALUES ('",line$UUID,"', 0, ",
                                        line$FOLD_CHANGE, ", ",
                                        line$MEAN_CORR, ", ",
                                        line$MEDIAN_CORR,");"))
  discard <- dbExecute(conn_new, paste0("INSERT INTO characteristics(",
                                        "UUID, PARTIAL, FOLD_CHANGE, ",
                                        "MEAN_CORR, MEDIAN_CORR) ",
                                        "VALUES ('",line$UUID,"', 1, ",
                                        line$FOLD_CHANGE_PARTIAL, ", ",
                                        line$MEAN_CORR_PARTIAL, ", ",
                                        line$MEDIAN_CORR_PARTIAL,");"))
}
for(i in 1:nrow(results_partial)) {
  line <- results_partial[i,,drop=FALSE]
  discard <- dbExecute(conn_new, paste0("INSERT INTO characteristics",
                                        "(UUID, PARTIAL) ",
                                        "VALUES ('",line$UUID,"', 0);"))
  discard <- dbExecute(conn_new, paste0("INSERT INTO characteristics",
                                        "(UUID, PARTIAL) ",
                                        "VALUES ('",line$UUID,"', 1);"))
}

# dbGetQuery(conn, "SELECT * FROM datasets LIMIT 10;")
# dbGetQuery(conn, "SELECT * FROM characteristics LIMIT 10;")
# dbGetQuery(conn, "SELECT * FROM results LIMIT 10;")

dbDisconnect(conn)
dbDisconnect(conn_new)
