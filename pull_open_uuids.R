source("path_fix.R")

library(RSQLite)
library(dplyr)
library(optparse)

option_list = list(
  make_option(c("--p"), type = "numeric", default = NULL,
              help = "number of features", metavar = "numeric"),
  make_option(c("--method"), type = "character", default = NULL,
              help = "differential abundance calling method", metavar = "character"),
  make_option(c("--corrp"), type = "numeric", default = 0,
              help = "correlation of features flag", metavar = "numeric"),
  make_option(c("--file"), type = "character", default = NULL,
              help = "save file name", metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

p <- opt$p
method <- opt$method
corrp <- opt$corrp
file <- opt$file

if(is.null(file) | file == "") {
  stop("Invalid save file name!\n")
}

# Testing
# p <- 1000
# method <- "DESeq2"
# corrp <- 0
# file <- "temp.txt"

if(!(method %in% c("ALDEx2", "DESeq2", "MAST", "NBGLM", "scran"))) {
  stop("Unrecognized method:", method, "\n")
}

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

all_uuids <- dbGetQuery(conn, paste0("SELECT UUID FROM datasets WHERE P = ",p," AND CORRP = ",corrp,";"))$UUID
wishlist <- data.frame(uuid = c(), baseline = c(), partial_info = c(), result_type = c())
for(u in 1:length(all_uuids)) {
  # Check for each of 4 fits
  cat(paste0("UUID ", u, " / ", length(all_uuids), "\n"))
  uuid <- all_uuids[u]
  res <- dbGetQuery(conn, paste0("SELECT * FROM results ",
                                 "WHERE UUID = '",uuid,"' AND ",
                                 "METHOD = '",method,"' AND ",
                                 "RESULT IS NOT NULL;"))
  # There are 8 combos of fit we want for each method
  #   1) 'self' or 'threshold' as differential abundance reference
  #   2) partial info = 0 or 1
  #   3) result type = tpr or fpr
  #
  # Look for any missing combos and add them to the wishlist
  for(ref in c("self", "threshold")) {
    for(par in c(0, 1)) {
      for(rt in c("tpr", "fpr")) {
        if(res %>% filter(BASELINE == ref & PARTIAL_INFO == par & RESULT_TYPE == rt) %>% count() %>% pull(n) == 0) {
          wishlist <- rbind(wishlist,
                            data.frame(uuid = uuid, baseline = ref, partial_info = par, result_type = rt))
        }
      }
    }
  }
}

dbDisconnect(conn)

wishlist <- wishlist %>%
  select(-result_type) %>%
  distinct()

write.table(wishlist, file = paste0(file, "_", p, "_", corrp, "_", method, ".txt"))
