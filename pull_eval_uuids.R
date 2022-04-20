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
file <- paste0("input_eval_", p, ".txt")

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

  # Look for any missing combos and add them to the wishlist
  # This is pretty run-specific!
  for(ref in c("oracle")) {
    for(par in c(0)) {
      #for(method in c("ALDEx2", "DESeq2", "scran", "edgeR", "edgeR_TMM")) {
      for(method in c("ANCOMBC")) {
        #for(obs in c("cpm")) {
        for(obs in c("relative_abundances")) {
          if(res %>% filter(BASELINE_TYPE == ref & PARTIAL_INFO == par & METHOD == method & OBSERVED_TYPE == obs) %>% count() %>% pull(n) == 0) {
            wishlist <- rbind(wishlist,
                              data.frame(ID = counter,
                                         uuid = uuid,
                                         baseline = ref,
                                         partial_info = par,
                                         method = method,
                                         observed_type = obs))
            counter <- counter + 1
          }
        }
      }
    }
  }
}

dbDisconnect(conn)

write.table(wishlist, file)
