source("path_fix.R")

library(RSQLite)
library(codaDE)
library(tidyverse)

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

# What I just updated
res1 <- dbGetQuery(conn, "SELECT UUID, METHOD, PARTIAL_INFO, BASELINE_TYPE, OBSERVED_TYPE, FDR, BETA FROM results WHERE FDR = 0.05")

# Everything else
res2 <- dbGetQuery(conn, "SELECT UUID, METHOD, PARTIAL_INFO, BASELINE_TYPE, OBSERVED_TYPE, FDR, BETA FROM results WHERE FDR != 0.05")

res <- res1 %>%
  left_join(res2, by = c("UUID", "METHOD", "PARTIAL_INFO", "BASELINE_TYPE", "OBSERVED_TYPE"))

# Remove these - they've been updated but the placeholders remain
remove <- res %>%
  filter(is.na(res$FDR.y))

ops <- 0
for(i in 1:nrow(remove)) {
  if(i %% 100 == 0) {
    cat(paste0("Iteration ", i, " / ", nrow(remove), "\n"))
  }
  ops <- ops + dbExecute(conn, paste0("DELETE FROM results WHERE UUID='", remove$UUID[i], "'",
                                                          " AND METHOD='", remove$METHOD[i], "'",
                                                          " AND PARTIAL_INFO=0",
                                                          " AND BASELINE_TYPE='", remove$BASELINE_TYPE[i], "'",
                                                          " AND OBSERVED_TYPE='", remove$OBSERVED_TYPE[i], "'",
                                                          " AND FDR IS NULL"))
}

dbDisconnect(conn)

cat(paste0(ops, " rows updated\n"))

