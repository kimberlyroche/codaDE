ranked_features <- NULL
for(result_type in c("TPR", "FPR")) {
  for(method in c("ALDEx2", "DESeq2", "edgeR_TMM", "scran")) {
    top_features <- readRDS(file.path("output",
                                      "predictive_fits",
                                      "oracle",
                                      "regression_plain",
                                      paste0("oracle_",result_type,"_",method,"_variable-importance.rds")))
    top_features$feature <- rownames(top_features)
    rownames(top_features) <- NULL
    ranked_features <- rbind(ranked_features,
                             cbind(data.frame(result_type = result_type,
                                              method = method,
                                              rank = 1:5),
                                   top_features %>%
                                     arrange(desc(Overall)) %>%
                                     slice(1:5) %>%
                                     mutate(importance = round(Overall/max(Overall), 2)) %>%
                                     select(feature, importance)))
  }
}

# Top features associated with TPR
ranked_features_tpr <- ranked_features %>%
  filter(result_type == "TPR") %>%
  arrange(method, rank)

ranked_features_fpr <- ranked_features %>%
  filter(result_type == "FPR") %>%
  arrange(method, rank)

ranked_features_tpr

ranked_features_fpr

# Calculate the Spearman's correlation between top features and outcomes
conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
query <- paste0("SELECT DISTINCT partB.UUID, partA.METHOD, partA.TPR, partA.FPR, ",
                paste0(paste0("partB.", unique(c(ranked_features_tpr$feature, ranked_features_fpr$feature))), collapse = ", "),
                " FROM (SELECT UUID, METHOD, TPR, FPR FROM results WHERE METHOD IN ('ALDEx2', 'DESeq2', 'edgeR_TMM', 'scran') AND
PARTIAL_INFO=0 AND OBSERVED_TYPE='relative_abundances' AND BASELINE_TYPE='oracle') AS partA
LEFT JOIN (SELECT * FROM characteristics WHERE TYPE='relative_abundances' AND PARTIAL=0) AS partB ON partA.UUID=partB.UUID")
res <- dbGetQuery(conn, query)
dbDisconnect(conn)

res <- res %>%
  filter(complete.cases(.))

ranked_features_tpr$rho <- NA
for(f in unique(ranked_features_tpr$feature)) {
  rho <- cor(res[,colnames(res) == f], res$TPR, method = "spearman")
  ranked_features_tpr$rho[ranked_features_tpr$feature == f] <- rho
}

ranked_features_fpr$rho <- NA
for(f in unique(ranked_features_fpr$feature)) {
  rho <- cor(res[,colnames(res) == f], 1-res$FPR, method = "spearman")
  ranked_features_fpr$rho[ranked_features_fpr$feature == f] <- rho
}

ranked_features_tpr

ranked_features_fpr

