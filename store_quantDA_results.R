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
file_list <- list.files(dir, pattern = "output_quantDA")
for(file in file_list) {
  results <- read.table(file.path(dir, file), header = TRUE)
  for(i in 1:nrow(results)) {
    job <- results[i,]
    updates <- updates + dbExecute(conn,
                                   paste0("INSERT OR REPLACE INTO da_realized(UUID, MEASURED_BY, CALLS) ",
                                          "VALUES('", job$uuid, "', '",job$measured_by,"', '", job$calls, "');"))
  }
  cat(paste0("Succeeded on ", updates, " rows\n"))
}

dbDisconnect(conn)
