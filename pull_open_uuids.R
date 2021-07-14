source("path_fix.R")

library(RSQLite)
library(dplyr)
library(optparse)

option_list = list(
  make_option(c("--p"), type = "numeric", default = NULL,
              help = "number of features", metavar = "numeric"),
  make_option(c("--file"), type = "character", default = NULL,
              help = "save file name", metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

p <- opt$p
file <- opt$file

if(is.null(file) | file == "") {
  stop("Invalid save file name!\n")
}

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

all_uuids <- dbGetQuery(conn, paste0("SELECT UUID FROM datasets WHERE P = ",p,";"))$UUID
wishlist <- data.frame(uuid = c(), baseline = c(), partial_info = c(), result_type = c())
counter <- 1
for(u in 1:length(all_uuids)) {
  if(u %% 100 == 0) {
    cat(paste0("UUID ", u, " / ", length(all_uuids), "\n"))
  }
  uuid <- all_uuids[u]
  res <- dbGetQuery(conn, paste0("SELECT * FROM results ",
                                 "WHERE UUID = '",uuid,"' AND ",
                                 "CALLS IS NOT NULL;"))
  # There are 16 combos of fit we want for each dataset
  #   1) 'self' or 'threshold' as differential abundance reference
  #   2) partial info = 0 or 1
  #   3) method = ALDEx2, DESeq2, MAST, scran
  #
  # Look for any missing combos and add them to the wishlist
  for(ref in c("self", "threshold")) {
    for(par in c(0, 1)) {
      for(method in c("ALDEx2", "DESeq2", "MAST", "scran")) {
        if(res %>% filter(BASELINE_TYPE == ref & PARTIAL_INFO == par & METHOD == method) %>% count() %>% pull(n) == 0) {
          wishlist <- rbind(wishlist,
                            data.frame(ID = counter, uuid = uuid, baseline = ref, partial_info = par, method = method))
          counter <- counter + 1
        }
      }
    }
  }
}

dbDisconnect(conn)

write.table(wishlist, file = paste0(file, "_", p, ".txt"))
