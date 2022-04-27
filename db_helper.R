source("path_fix.R")

library(RSQLite)
library(codaDE)
library(tidyverse)

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

# res <- dbGetQuery(conn, "SELECT * FROM results_backup WHERE METHOD='DESeq2_control' AND BASELINE_TYPE='oracle' AND OBSERVED_TYPE='relative_abundances' AND PARTIAL_INFO=0")
res <- dbGetQuery(conn, "SELECT * FROM results_backup WHERE METHOD='edgeR' AND BASELINE_TYPE='oracle' AND OBSERVED_TYPE='relative_abundances' AND PARTIAL_INFO=0")

ins <- 0
for(i in 1:nrow(res)) {
  if(is.na(res$TPR[i]) | is.na(res$FPR[i])) next;
  ins <- ins + dbExecute(conn, paste0("INSERT OR IGNORE INTO results(UUID, METHOD, PARTIAL_INFO, BASELINE_TYPE, OBSERVED_TYPE, BASELINE_CALLS, CALLS, BASELINE_BETAS, BETAS, TPR, FPR, FDR, BETA) VALUES(",
    "'", res$UUID[i], "', ",
    "'", res$METHOD[i], "', ",
    res$PARTIAL_INFO[i], ", ",
    "'", res$BASELINE_TYPE[i], "', ",
    "'", res$OBSERVED_TYPE[i], "', ",
    ifelse(is.na(res$BASELINE_CALLS[i]), "NULL", paste0("'", res$BASELINE_CALLS[i], "'")), ", ",
    "'", res$CALLS[i], "', ",
    "NULL, ",
    "NULL, ",
    res$TPR[i], ", ",
    res$FPR[i], ", ",
    "0.05, ",
    "-1);"))
}
cat(paste0("Inserted ", ins, " rows\n"))

dbDisconnect(conn)

