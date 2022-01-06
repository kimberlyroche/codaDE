source("path_fix.R")

library(codaDE)
library(RSQLite)
library(dplyr)
library(ggplot2)

threshold <- 1.33
percent_agreement <- NULL

for(this_method in c("ALDEx2", "DESeq2", "scran")) {
  for(this_P in c(100, 1000, 5000)) {
    cat(paste0("Method: ", this_method, ", P: ", this_P, "\n"))
    conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
    # This takes upwards of 30 seconds for one METHOD-P combo
    results <- dbGetQuery(conn, paste0("SELECT datasets.UUID, METHOD, ",
                                       "results.BASELINE_CALLS, P, PARTIAL_INFO, ",
                                       "BASELINE_TYPE FROM results LEFT JOIN datasets ",
                                       "ON results.UUID=datasets.UUID WHERE ",
                                       "METHOD='", this_method, "' AND P=",this_P," AND ",
                                       "PARTIAL_INFO=0 AND BASELINE_TYPE='self'"))
    results_fc <- dbGetQuery(conn, paste0("SELECT UUID, CALLS FROM da_realized WHERE UUID IN ('",
                                          paste0(results$UUID, collapse="', '"),
                                          "') AND MEASURED_BY='fc_",threshold,"'"))
    
    joined_results <- results %>%
      inner_join(results_fc, by = "UUID")
    
    for(i in 1:nrow(joined_results)) {
      x <- as.numeric(strsplit(joined_results[i,]$BASELINE_CALLS, ";")[[1]])
      y <- as.numeric(strsplit(joined_results[i,]$CALLS, ";")[[1]])
      z <- rep(1, length(x))
      z[x < 0.05] <- 0
      percent_agreement <- rbind(percent_agreement,
                                 data.frame(method = this_method,
                                            P = this_P,
                                            percent_agree = sum(y == z)/length(y)))
    }
  }
}

p <- ggplot(percent_agreement, aes(x = percent_agree)) +
  geom_histogram(color = "white") +
  facet_wrap(. ~ method + P) +
  theme_bw() +
  labs(x = paste0("percent agreement: ", threshold, "x FC and calls on absolute counts"))
ggsave(file.path("output", "images", paste0("FC_vs_baseline_calls_", threshold, ".png")),
       p,
       dpi = 100,
       units = "in",
       height = 9,
       width = 9)

