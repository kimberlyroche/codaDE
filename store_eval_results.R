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

insertions <- 0
file_list <- list.files(dir, pattern = "output_eval")
for(file in file_list) {
  results <- read.table(file.path(dir, file), sep = "\t", header = TRUE)
  # Fix messed up headers
  colnames(results) <- c("id", "uuid", "baseline", "partial_info", "method", "observed_type", "baseline_calls", "baseline_betas", "calls", "betas")
  for(i in 1:nrow(results)) {
    job <- results[i,]
    if(!any(is.na(job %>% select(-c(baseline_calls, baseline_betas, betas))))) {
      insertions <- insertions + dbExecute(conn, paste0("INSERT OR IGNORE INTO RESULTS(UUID, METHOD, PARTIAL_INFO, BASELINE_TYPE, OBSERVED_TYPE, ",
                                                        "BASELINE_CALLS, BASELINE_BETAS, CALLS, BETAS) VALUES (",
                                                        "'",job$uuid,"', ",
                                                        "'",job$method,"', ",
                                                        job$partial_info,", ",
                                                        "'",job$baseline,"', ",
                                                        "'",job$observed_type,"', ",
                                                        ifelse(is.na(job$baseline_calls), "NULL, ", paste0("'", job$baseline_calls,"', ")),
                                                        ifelse(is.na(job$baseline_betas), "NULL, ", paste0("'", job$baseline_betas,"', ")),
                                                        "'",job$calls,"', ",
                                                        ifelse(is.na(job$betas), "NULL, ", paste0("'",job$betas,"')"))))
    }
  }
  cat(paste0("Succeeded on ", insertions, " rows\n"))
}

dbDisconnect(conn)
