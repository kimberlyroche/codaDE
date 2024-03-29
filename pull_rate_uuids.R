source("path_fix.R")

library(RSQLite)
library(dplyr)
library(tidyr)
library(optparse)

option_list = list(
  make_option(c("--p"), type = "numeric", default = NULL,
              help = "number of features", metavar = "numeric"),
  make_option(c("--alpha"), type = "numeric", default = 0.05,
              help = "FDR alpha", metavar = "numeric"),
  make_option(c("--beta"), type = "numeric", default = NULL,
              help = "effect size threshold", metavar = "numeric")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

p <- opt$p
alpha <- opt$alpha
beta <- opt$beta
file <- paste0("input_rate_", p, ".txt")

#methods_list <- c("ALDEx2", "ANCOMBC", "DESeq2", "scran", "edgeR", "edgeR_TMM", "DESeq2_control")
methods_list <- c("scran")

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

all_uuids <- dbGetQuery(conn,
                        paste0("SELECT UUID FROM datasets WHERE P = ",p,";"))$UUID
wishlist <- data.frame(ID = c(),
                       uuid = c(),
                       baseline = c(),
                       partial_info = c(),
                       method = c(),
                       type = c(),
                       baseline_calls = c(),
                       calls = c())
counter <- 1
for(u in 1:length(all_uuids)) {
  if(u %% 100 == 0) {
    cat(paste0("UUID ", u, " / ", length(all_uuids), "\n"))
  }
  uuid <- all_uuids[u]
  # Inclusion scenario 1: missing TPR or FPR for existing rows
  # Pull all recorded results
  res <- dbGetQuery(conn, paste0("SELECT results.UUID, METHOD, PARTIAL_INFO, ",
                                 "BASELINE_TYPE, OBSERVED_TYPE, ",
                                 "datasets.BASELINE_CALLS AS ORACLE_CALLS, ",
                                 "datasets.BASELINE_BETAS AS ORACLE_BETAS, ",
                                 "results.BASELINE_CALLS, results.BASELINE_BETAS, ",
                                 "CALLS, BETAS, TPR, FPR, FDR, BETA ",
                                 "FROM results LEFT JOIN datasets ",
                                 "ON results.UUID = datasets.UUID ",
                                 "WHERE results.UUID = '",uuid,"' AND ",
                                 "CALLS IS NOT NULL;"))

  # Retain those where TPR or FPR have not been recorded
  res2 <- res %>%
    filter((is.na(TPR) | is.na(FPR)) & FDR==alpha)

  # Inclusion scenario 2: rows for this FDR associated with existing results
  # res3 <- res %>%
  #   filter(FDR != alpha) %>%
  #   mutate(FDR = alpha)
  res3 <- res %>%
    filter(is.na(FDR) | FDR != alpha) %>%
    left_join(res %>%
                filter(FDR == alpha) %>%
                select(-c(ORACLE_CALLS, ORACLE_BETAS, BASELINE_CALLS, BASELINE_BETAS, CALLS, BETAS)),
              by = c("UUID" = "UUID",
                     "METHOD" = "METHOD",
                     "PARTIAL_INFO" = "PARTIAL_INFO",
                     "BASELINE_TYPE" = "BASELINE_TYPE",
                     "OBSERVED_TYPE" = "OBSERVED_TYPE")) %>%
    filter(is.na(TPR.y) | is.na(FPR.y) | is.na(FDR.y)) %>%
    mutate(FDR = alpha) %>%
    select(-c(FDR.x, TPR.y, FPR.y, FDR.y, BETA.x, BETA.y))
  
  res <- rbind(res2, res3)

  if(nrow(res) > 0) {

    if(alpha != 0.05 & !is.null(beta)) {
      res$BETA <- beta
    } else {
      res$BETA <- -1
    }
    
    if(nrow(res) > 0) {
      for(i in 1:nrow(res)) {
        job <- res[i,]
        if(job$METHOD %in% methods_list & job$BASELINE_TYPE == "oracle" & job$OBSERVED_TYPE == "relative_abundances") {
          wishlist <- rbind(wishlist,
                            data.frame(ID = counter,
                                       uuid = job$UUID,
                                       baseline = job$BASELINE_TYPE,
                                       partial_info = job$PARTIAL_INFO,
                                       method = job$METHOD,
                                       type = job$OBSERVED_TYPE,
                                       baseline_calls = ifelse(job$BASELINE_TYPE == "oracle", job$ORACLE_CALLS, job$BASELINE_CALLS),
                                       baseline_betas = ifelse(job$BASELINE_TYPE == "oracle", job$ORACLE_BETAS, job$BASELINE_BETAS),
                                       calls = job$CALLS,
                                       betas = job$BETAS,
                                       fdr = job$FDR,
                                       beta = job$BETA))
          counter <- counter + 1
        }
      }
    }
  }
}

dbDisconnect(conn)

write.table(wishlist, file)
