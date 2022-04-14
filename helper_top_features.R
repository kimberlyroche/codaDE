ranked_features <- NULL

for(result_type in c("TPR", "FPR")) {
  for(method in c("ALDEx2", "DESeq2", "scran")) {
    top_features <- readRDS(file.path("output",
                                      "predictive_fits",
                                      "oracle",
                                      "regression_cpm",
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

ranked_features %>%
  filter(result_type == "TPR") %>%
  arrange(method, rank)

ranked_features %>%
  filter(result_type == "FPR") %>%
  arrange(method, rank)
