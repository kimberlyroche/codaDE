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
  results <- read.table(file.path(dir, file), header = TRUE)
  for(i in 1:nrow(results)) {
    job <- results[i,]
    if(!any(is.na(job))) {
      insertions <- insertions + dbExecute(conn, paste0(
        "INSERT INTO datasets(UUID, P, CORRP, LOG_MEAN, ",
                             "PERTURBATION, REP_NOISE, ",
                             "FC_ABSOLUTE, FC_RELATIVE, ",
                             "FC_PARTIAL, BASELINE_CALLS) ",
        "VALUES(",
        "'", job$UUID, "', ",
        job$P, ", ",
        job$CORRP, ", ",
        job$LOG_MEAN, ", ",
        job$PERTURBATION, ", ",
        job$REP_NOISE, ", ",
        job$FC_ABSOLUTE, ", ",
        job$FC_RELATIVE, ", ",
        job$FC_PARTIAL, ", ",
        "'", job$BASELINE_CALLS, "'",
        ")"
      ))
    }
  }
  cat(paste0("Succeeded on ", insertions, " rows\n"))
}

dbDisconnect(conn)
