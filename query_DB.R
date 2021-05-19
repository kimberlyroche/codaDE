source("path_fix.R")

library(RSQLite)
library(optparse)

option_list = list(
  make_option(c("--query"), type = "character", default = "results",
              help = "database query", metavar = "numeric")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
dbGetQuery(conn, opt$query)
dbDisconnect(conn)
