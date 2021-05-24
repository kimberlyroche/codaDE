source("path_fix.R")

library(codaDE)
library(RSQLite)

# Create DB (if doesn't exist?)
conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

discard <- dbExecute(conn, paste0("CREATE TABLE datasets(",
                                  "UUID VARCHAR(36) PRIMARY KEY,",
                                  "P INTEGER,",
                                  "CORRP INTEGER,",
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
                                  "PARTIAL_INFO INT,",
                                  "RESULT REAL,",
                                  "RESULT_TYPE VARCHAR(24),",
                                  "FIT_IN_PROGRESS INTEGER);"))

dbDisconnect(conn)
