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
file <- paste0("input_char_", p, ".txt")

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

all_uuids <- dbGetQuery(conn, paste0("SELECT UUID FROM datasets WHERE P = ",p,";"))$UUID
wishlist <- data.frame(uuid = c(), baseline = c(), partial_info = c(), result_type = c())
counter <- 1
for(u in 1:length(all_uuids)) {
  if(u %% 100 == 0) {
    cat(paste0("UUID ", u, " / ", length(all_uuids), "\n"))
  }
  uuid <- all_uuids[u]
  for(par in c(0, 1)) {
    res <- dbGetQuery(conn, paste0("SELECT * FROM characteristics ",
                                   "WHERE UUID = '",uuid,"' AND ",
                                   "PARTIAL=", par, ";"))
    # Look for any missing fields
    if(nrow(res) == 0 || any(is.na(res %>% select(-c(UUID, PARTIAL))))) {
      wishlist <- rbind(wishlist,
                        data.frame(ID = counter, uuid = uuid, partial_info = par))
      counter <- counter + 1
    }
  }
}

dbDisconnect(conn)

write.table(wishlist, file)
