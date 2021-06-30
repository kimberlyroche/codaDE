# This helper script will cycle through all datasets and update the PERCENT_DIFF
# field in the `datasets` table in the DB.

source("path_fix.R")

library(RSQLite)
library(ggplot2)

library(optparse)

option_list = list(
  make_option(c("--p"), type = "numeric", default = 100,
              help = "number of features", metavar = "numeric"),
  make_option(c("--corrp"), type = "numeric", default = 0,
              help = "correlation flag", metavar = "numeric")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

p <- opt$p
corrp <- opt$corrp

# Testing
# p <- 100
# corrp <- 0

# Pull UUIDs matching P, CORRP characterstics
conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
res <- dbGetQuery(conn, paste0("SELECT UUID FROM datasets WHERE P=",p," AND CORRP=",corrp))
updates <- 0
for(i in 1:length(res$UUID)) {
  cat(paste0("Parsing dataset ", i, " / ", length(res$UUID), "\n"))
  uuid <- res$UUID[i]
  data <- readRDS(file.path("output", "datasets", paste0(uuid, ".rds")))
  # Estimate % DE features
  M <- data$simulation$abundances
  prop_de <- sum(calc_threshold_DA(M) == 0) / ncol(M)
  discard <- dbExecute(conn, paste0("UPDATE datasets SET PERCENT_DIFF=",prop_de," WHERE UUID='",uuid,"'"))
  updates <- updates + discard
}
dbDisconnect(conn)

cat(paste0("Updated ", updates, " datasets\n"))
