source("path_fix.R")

library(RSQLite)
library(codaDE)
library(tidyverse)
library(optparse)

# -----------------------------------------------------------------------------
#   PART 1
# -----------------------------------------------------------------------------

if(FALSE) {
  option_list = list(
    make_option(c("--start"), type = "numeric", default = 1,
                help = "start row index", metavar = "numeric"),
    make_option(c("--end"), type = "numeric", default = Inf,
                help = "end row index", metavar = "numeric")
  );

  opt_parser = OptionParser(option_list = option_list);
  opt = parse_args(opt_parser);

  start <- opt$start
  end <- opt$end

  output_fn <- file.path("temp", paste0("betas_", start, "-", end, ".txt"))

  # Initialize output
  output_file <- file(output_fn)
  writeLines(paste0(c("uuid",
                      "betas"), collapse = "\t"),
             output_file)
  close(output_file)

  # Pull all files in output/datasets
  files <- list.files(file.path("output", "datasets"))
  for(f in start:end) {
    cat(paste0("File ", f, " / ", length(files), "\n"))
    file <- files[f]
    sim_data <- readRDS(file.path("output", "datasets", file))$simulation
    calls_obj <- call_DA_NB(sim_data$abundances, sim_data$groups)
    uuid <- str_replace(file, "\\.rds", "")
    betas <- paste(calls_obj$beta, collapse = ";")
    write_delim(data.frame(uuid = uuid, betas = betas), output_fn, delim = "\t", append = TRUE)
  }
}

# -----------------------------------------------------------------------------
#   PART 2
# -----------------------------------------------------------------------------

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

# Pull all files in output/datasets
added <- 0
files <- list.files("temp")
for(f in 1:length(files)) {
  file <- files[f]
  results <- read.delim(file.path("temp", file), sep = "\t", header = TRUE)
  for(i in 1:nrow(results)) {
    res <- results[i,]
    query <- dbExecute(conn, paste0("UPDATE datasets SET BASELINE_BETAS='", paste(res$beta, collapse = ";"), "' WHERE UUID='", res$uuid,"'"))
    added <- added + query
    if(added %% 100 == 0) {
      cat(paste0("Added row ", added, "\n"))
    }
  }
}

dbDisconnect(conn)
