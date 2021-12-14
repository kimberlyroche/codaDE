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

measured_by <- "fc_1.5"
p <- opt$p
file <- paste0("input_quantDA_", p, ".txt")

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

all_uuids <- dbGetQuery(conn, paste0("SELECT UUID FROM datasets WHERE P = ",p,";"))$UUID
wishlist <- data.frame(uuid = c())
counter <- 1
for(u in 1:length(all_uuids)) {
  if(u %% 100 == 0) {
    cat(paste0("UUID ", u, " / ", length(all_uuids), "\n"))
  }
  uuid <- all_uuids[u]
    res <- dbGetQuery(conn, paste0("SELECT * FROM da_realized ",
                                   "WHERE UUID = '",uuid,"' AND ",
                                   "MEASURED_BY='", measured_by, "';"))
    # Look for any missing fields
    if(nrow(res) == 0 || any(is.na(res %>% select(PERCENT_DIFF_REALIZ)))) {
      wishlist <- rbind(wishlist,
                        data.frame(ID = counter, uuid = uuid))
      counter <- counter + 1
    }
}

dbDisconnect(conn)

write.table(wishlist, file)
