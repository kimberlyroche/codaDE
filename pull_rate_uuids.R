source("path_fix.R")

library(RSQLite)
library(dplyr)
library(optparse)

option_list = list(
  make_option(c("--p"), type = "numeric", default = NULL,
              help = "number of features", metavar = "numeric")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

p <- opt$p
file <- paste0("input_rate_", p, ".txt")

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

all_uuids <- dbGetQuery(conn,
                        paste0("SELECT UUID FROM datasets WHERE P = ",p,";"))$UUID
wishlist <- data.frame(ID = c(),
                       uuid = c(),
                       baseline = c(),
                       partial_info = c(),
                       method = c(),
                       baseline_calls = c(),
                       calls = c())
counter <- 1
for(u in 1:length(all_uuids)) {
  if(u %% 100 == 0) {
    cat(paste0("UUID ", u, " / ", length(all_uuids), "\n"))
  }
  uuid <- all_uuids[u]
  res <- dbGetQuery(conn, paste0("SELECT results.UUID, METHOD, PARTIAL_INFO, ",
                                 "BASELINE_TYPE, ",
                                 "datasets.BASELINE_CALLS AS ORACLE_CALLS, ",
                                 "results.BASELINE_CALLS, CALLS, TPR, FPR ",
                                 "FROM results LEFT JOIN datasets ",
                                 "ON results.UUID = datasets.UUID ",
                                 "WHERE results.UUID = '",uuid,"' AND ",
                                 "CALLS IS NOT NULL AND ",
                                 "(TPR IS NULL OR FPR IS NULL);"))
  if(nrow(res) > 0) {
    for(i in 1:nrow(res)) {
      job <- res[i,]
      wishlist <- rbind(wishlist,
                        data.frame(ID = counter,
                                   uuid = job$UUID,
                                   baseline = job$BASELINE_TYPE,
                                   partial_info = job$PARTIAL_INFO,
                                   method = job$METHOD,
                                   baseline_calls = ifelse(job$BASELINE_TYPE == "oracle", job$ORACLE_CALLS, job$BASELINE_CALLS),
                                   calls = job$CALLS))
      counter <- counter + 1
    }
  }
}

dbDisconnect(conn)

write.table(wishlist, file)
