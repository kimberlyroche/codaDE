source("path_fix.R")

library(RSQLite)
library(optparse)

option_list = list(
  make_option(c("--table"), type = "character", default = "results",
              help = "database table name", metavar = "numeric")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

table <- opt$table
if(!(table %in% c("datasets", "results"))) {
  stop(paste0("Invalid table name: ", table, "!\n"))
}

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
dbGetQuery(conn, paste0("SELECT * FROM ",table,";"))
dbDisconnect(conn)
