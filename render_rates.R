source("path_fix.R")

library(codaDE)
library(RSQLite)

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

ps <- c(100, 1000, 5000)
partials <- c(0, 1)
references <- c("oracle", "self")

for(p in ps) {
  for(partial in partials) {
    for(reference in references) {
      cat(paste0("Evaluating P=", p, ", PARTIAL=", partial, ", REF=", reference, "\n"))

      res <- dbGetQuery(conn, paste0("SELECT ",
                                     "datasets.UUID AS UUID, ",
                                     "results.METHOD, results.PARTIAL_INFO, ",
                                     "BASELINE_TYPE, ",
                                     "results.BASELINE_CALLS AS SELF_BASELINE, ",
                                     "datasets.BASELINE_CALLS AS ORACLE_BASELINE, ",
                                     "CALLS, TPR, FPR, PERCENT_DIFF ",
                                     "FROM results LEFT JOIN datasets ON ",
                                     "results.UUID=datasets.UUID WHERE ",
                                     "P=", p,
                                     " AND PARTIAL_INFO=", partial,
                                     " AND BASELINE_TYPE='",reference,"'",
                                     " AND CALLS IS NOT NULL"))
      
      # Manually calculate TPR and FPR
      missing_rows <- unique(c(which(is.na(res$TPR)),
                               which(is.na(res$FPR)),
                               which(is.na(res$PERCENT_DIFF))))
      for(mr_idx in 1:length(missing_rows)) {
        cat(paste0("Updating incomplete results ", mr_idx, " / ", length(missing_rows), "\n"))
        i <- missing_rows[mr_idx]
        if(res[i,]$BASELINE_TYPE == "oracle") {
          rates <- calc_DA_discrepancy(res[i,]$CALLS, res[i,]$ORACLE_BASELINE)
        } else {
          rates <- calc_DA_discrepancy(res[i,]$CALLS, res[i,]$SELF_BASELINE)
        }
        if(!is.na(rates$TPR) && !is.na(rates$FPR)) {
          discard <- dbExecute(conn, paste0("UPDATE results SET ",
                                            "TPR=", rates$TPR, ", ",
                                            "FPR=", rates$FPR,
                                            " WHERE ",
                                            "UUID='", res$UUID[i], "' AND ",
                                            "METHOD='", res$METHOD[i], "' AND ",
                                            "PARTIAL_INFO=", res$PARTIAL_INFO[i], " AND ",
                                            "BASELINE_TYPE='", res$BASELINE_TYPE[i], "';"))
        } else {
          cat("\tNo detectable differential abundance!\n")
        }
        discard <- dbExecute(conn, paste0("UPDATE datasets SET ",
                                          "PERCENT_DIFF=", rates$percent_DA,
                                          " WHERE UUID='", res$UUID[i], "';"))
      }
    }
  }
}
