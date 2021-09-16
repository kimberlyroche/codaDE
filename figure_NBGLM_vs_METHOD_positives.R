library(codaDE)
library(ggplot)
library(RSQLite)

dir.create(file.path("output", "images"), showWarnings = FALSE)

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
                               "AND BASELINE_TYPE='self' ",
                               "AND FC_ABSOLUTE <= 10 ",
                               "AND FC_ABSOLUTE >= 0.1;"))
      
# Strip "result-less" entries
res <- res %>%
  filter(!is.na(TPR) & !is.na(FPR))

dbDisconnect(conn)

results <- NULL

for(method in c("ALDEx2", "DESeq2", "scran")) {
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

pl <- ggplot(results, aes(x = oracle_total, y = self_total, fill = FPR)) +
  geom_point(size = 2, shape = 21) +
  facet_wrap(. ~ method, ncol = 4) +
  theme_bw() +
  scale_fill_gradient2(low = "#bbbbbb", high = "red") +
  labs(x = "Number of differential features (NB GLM)",
       y = "Number of differential features (this method)") +
  theme(axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)))
ggsave(file.path("output", "images", "methods_vs_NBGLM.png"),
       pl,
       dpi = 100,
       units = "in",
       height = 4,
       width = 11)

# ------------------------------------------------------------------------------
#   Print R-squared values
# ------------------------------------------------------------------------------

for(this_method in c("ALDEx2", "DESeq2", "scran")) {
  cat(paste0("R^2 ", this_method, ": ",
             round(cor(results %>% filter(method == this_method) %>% pull(oracle_total),
                       results %>% filter(method == this_method) %>% pull(self_total))**2, 3), "\n"))
}


