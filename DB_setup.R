source("path_fix.R")

library(codaDE)
library(RSQLite)

# Create DB (if doesn't exist?)
conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

discard <- dbExecute(conn, paste0("CREATE TABLE datasets(",
                                  "UUID VARCHAR(36) PRIMARY KEY,",
                                  "P INTEGER,",
                                  "CORRP REAL,",
                                  "TIMESTAMP VARCHAR(19),",
                                  "FOLD_CHANGE REAL,",
                                  "FOLD_CHANGE_PARTIAL REAL,",
                                  "MEAN_CORR REAL,",
                                  "MEAN_CORR_PARTIAL REAL,",
                                  "MEDIAN_CORR REAL,",
                                  "MEDIAN_CORR_PARTIAL REAL);"))
discard <- dbExecute(conn, paste0("CREATE TABLE results(",
                                  "RUN_ID INTEGER PRIMARY KEY AUTOINCREMENT,",
                                  "UUID VARCHAR(36),",
                                  "METHOD VARCHAR(64),",
                                  "PARTIAL_INFO INTEGER,",
                                  "RESULT REAL,",
                                  "RESULT_TYPE VARCHAR(24),",
                                  "FIT_IN_PROGRESS INTEGER);"))

# Add timestamp columns and update NAs
# discard <- dbExecute(conn, "ALTER TABLE datasets ADD COLUMN TIMESTAMP VARCHAR(19);")
# discard <- dbExecute(conn, paste0("UPDATE datasets SET TIMESTAMP = '", get_timestamp(), "' WHERE TIMESTAMP IS NULL;"))

# Add feature columns
# discard <- dbExecute(conn, "ALTER TABLE datasets ADD COLUMN FOLD_CHANGE REAL;")
# discard <- dbExecute(conn, "ALTER TABLE datasets ADD COLUMN FOLD_CHANGE_PARTIAL REAL;")
# discard <- dbExecute(conn, "ALTER TABLE datasets ADD COLUMN MEAN_CORR REAL;")
# discard <- dbExecute(conn, "ALTER TABLE datasets ADD COLUMN MEAN_CORR_PARTIAL REAL;")
# discard <- dbExecute(conn, "ALTER TABLE datasets ADD COLUMN MEDIAN_CORR REAL;")
# discard <- dbExecute(conn, "ALTER TABLE datasets ADD COLUMN MEDIAN_CORR_PARTIAL REAL;")

# Add/drop fit-in-progress column
# discard <- dbExecute(conn, "ALTER TABLE results ADD COLUMN FIT_IN_PROGRESS INT;")
# discard <- dbExecute(conn, "ALTER TABLE results DROP COLUMN FIT_IN_PROGRESS;")

# List tables
# dbListTables(conn)

# dbGetQuery(conn, "SELECT * FROM datasets;")
# dbGetQuery(conn, "SELECT * FROM results;")

# Flush DB
# discard <- dbExecute(conn, "DELETE FROM results WHERE TRUE;")

dbDisconnect(conn)
