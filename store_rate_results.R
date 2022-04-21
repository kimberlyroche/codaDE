source("path_fix.R")

library(RSQLite)
library(dplyr)
library(optparse)

option_list = list(
  make_option(c("--dir"), type = "character", default = "temp",
              help = "directory of evaluation output files", metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

dir <- opt$dir

if(!file.exists(dir)) {
  stop("No such output directory!\n")
}

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

updates <- 0
file_list <- list.files(dir, pattern = "output_rate")
for(file in file_list) {
  results <- read.table(file.path(dir, file), sep = "\t", header = TRUE)
  results <- results[complete.cases(results),]
  for(i in 1:nrow(results)) {
    job <- results[i,]
    this_round <- dbExecute(conn,
                            paste0("UPDATE RESULTS SET ",
                                   "TPR = ",job$tpr, ", ",
                                   "FPR = ", job$fpr, " ",
                                   "WHERE ",
                                   "UUID = '", job$uuid, "' AND ",
                                   "METHOD = '", job$method, "' AND ",
                                   "PARTIAL_INFO = ", job$partial_info, " AND ",
                                   "BASELINE_TYPE = '", job$baseline, "' AND ",
                                   "OBSERVED_TYPE = '", job$type, "' AND ",
                                   "FDR = ", job$fdr))
    if(this_round == 0) {
      # Attempt to insert new row
      # First we'll need to pull some extra info from comparable results
      res <- dbGetQuery(conn, paste0("SELECT BASELINE_CALLS, CALLS FROM results WHERE UUID='", job$uuid, "'",
                                     " AND METHOD='", job$method, "' AND PARTIAL_INFO=", job$partial_info,
                                     " AND BASELINE_TYPE='", job$baseline, "' AND OBSERVED_TYPE='", job$type, "' LIMIT 1"))
      this_round <- dbExecute(conn, paste0("INSERT OR IGNORE INTO RESULTS(UUID, METHOD, PARTIAL_INFO, BASELINE_TYPE, OBSERVED_TYPE, BASELINE_CALLS, CALLS, TPR, FPR, FDR) VALUES (",
                                                        "'",job$uuid,"', ",
                                                        "'",job$method,"', ",
                                                        job$partial_info,", ",
                                                        "'",job$baseline,"', ",
                                                        "'",job$type,"', ",
                                                        ifelse(is.na(res$BASELINE_CALLS), "NULL, ", paste0("'", res$BASELINE_CALLS,"', ")),
                                                        "'",res$CALLS,"', ",
                                                        job$tpr,", ",
                                                        job$fpr,", ",
                                                        job$fdr,")"))
    }
    updates <- updates + this_round
  }
  cat(paste0("Succeeded on ", updates, " rows\n"))
}

dbDisconnect(conn)
