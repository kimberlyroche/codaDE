library(codaDE)
library(ggplot)
library(RSQLite)

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

res <- dbGetQuery(conn, paste0("SELECT ",
                               "datasets.UUID AS UUID, ",
                               "METHOD, ",
                               "PARTIAL_INFO, ",
                               "BASELINE_TYPE, ",
                               "datasets.BASELINE_CALLS AS ORACLE_BASELINE, ",
                               "CALLS, ",
                               "results.BASELINE_CALLS AS SELF_BASELINE, ",
                               "P, ",
                               "CORRP, ",
                               "LOG_MEAN, ",
                               "PERTURBATION, ",
                               "REP_NOISE, ",
                               "FC_ABSOLUTE, ",
                               "FC_RELATIVE, ",
                               "FC_PARTIAL, ",
                               "MED_ABS_TOTAL, ",
                               "MED_REL_TOTAL, ",
                               "PERCENT_DIFF_REALIZ, ",
                               "TPR, ",
                               "FPR ",
                               "FROM results LEFT JOIN datasets ON ",
                               "results.UUID=datasets.UUID ",
                               "WHERE PARTIAL_INFO=0 ",
                               # "AND P=100 ",
                               "AND BASELINE_TYPE='self' ",
                               "AND FC_ABSOLUTE <= 10 ",
                               "AND FC_ABSOLUTE >= 0.1;"))
      
# Strip "result-less" entries
res <- res %>%
  filter(!is.na(TPR) & !is.na(FPR))

dbDisconnect(conn)

results <- NULL

for(method in c("ALDEx2", "DESeq2", "MAST", "scran")) {
  cat(paste0("Making calls for method ", method, "\n"))
  # Filter to method of interest
  res_method <- res %>%
    filter(METHOD == method)
  # Pull oracle 'yes' total
  NB_calls <- unname(sapply(res_method$ORACLE_BASELINE,
                            function(x) {
                              sum(p.adjust(as.numeric(strsplit(x, ";")[[1]]), method = "BH") < 0.05)
                            }))
  # Pull method 'yes' total
  method_calls <- unname(sapply(res_method$SELF_BASELINE,
                                function(x) {
                                  sum(p.adjust(as.numeric(strsplit(x, ";")[[1]]), method = "BH") < 0.05)
                                }))
  results <- rbind(results,
                   data.frame(method = method,
                              oracle_total = NB_calls,
                              self_total = method_calls,
                              FPR = res_method$FPR))
}

test <- results[sample(1:nrow(results), size = 10000),]

ggplot(test, aes(x = oracle_total, y = self_total, fill = FPR)) +
  geom_point(size = 2, shape = 21) +
  facet_wrap(. ~ method, ncol = 4) +
  theme_bw() +
  scale_fill_gradient2(low = "#bbbbbb", high = "red") +
  labs(x = "% differential features (NB GLM)",
       y = "% differential features (this method)")







