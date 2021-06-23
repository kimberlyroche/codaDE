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
file_list <- list.files(dir)
for(file in file_list) {
  results <- read.table(file.path(dir, file))
  for(i in 1:nrow(results)) {
    job <- results[i,]
    insertions <- insertions + dbExecute(conn, paste0("INSERT INTO RESULTS(UUID, METHOD, PARTIAL_INFO, BASELINE, RESULT_TYPE, RESULT) VALUES (",
                                                      "'",job$uuid,"', ",
                                                      "'",job$method,"', ",
                                                      job$partial_info,", ",
                                                      "'",job$baseline,"', ",
                                                      "'tpr', ",
                                                      job$tpr,")"))
    insertions <- insertions + dbExecute(conn, paste0("INSERT INTO RESULTS(UUID, METHOD, PARTIAL_INFO, BASELINE, RESULT_TYPE, RESULT) VALUES (",
                                                      "'",job$uuid,"', ",
                                                      "'",job$method,"', ",
                                                      job$partial_info,", ",
                                                      "'",job$baseline,"', ",
                                                      "'fpr', ",
                                                      job$fpr,")"))
  }
  cat(paste0("Succeeded on ", insertions, " / ", nrow(results)*2, " rows\n"))
}

dbDisconnect(conn)
