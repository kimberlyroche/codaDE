library(tidyverse)

source("feature_map.R")
feature_map <- data.frame(symbol = names(feature_map),
                          print_name = unname(feature_map))

# Print all features
str_out <- ""
for(i in 1:nrow(feature_map)) {
  str_out <- paste0(str_out,
                    str_replace_all(feature_map[i,]$symbol, "_", "\\\\_"),
                    " & ",
                    feature_map[i,]$print_name,
                    " \\\\ \\hline \n")
}
writeLines(str_out, file.path("output", "scratch.txt"))

# Gather most important/informative features for prediction across methods
ranked_features <- NULL
for(result_type in c("TPR", "FPR")) {
  for(method in c("ALDEx2", "ANCOMBC", "DESeq2", "edgeR_TMM", "scran")) {
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
                                              rank = 1:3),
                                   top_features %>%
                                     arrange(desc(Overall)) %>%
                                     slice(1:3) %>%
                                     mutate(importance = round(Overall/max(Overall), 2)) %>%
                                     dplyr::select(feature, importance)))
  }
}

ranked_features$method <- factor(ranked_features$method)
levels(ranked_features$method)[2] <- "ANCOM-BC"
levels(ranked_features$method)[4] <- "edgeR (TMM)"

# Top features associated with TPR
ranked_features_tpr <- ranked_features %>%
  filter(result_type == "TPR") %>%
  arrange(method, rank)

ranked_features_fpr <- ranked_features %>%
  filter(result_type == "FPR") %>%
  arrange(method, rank)

temp <- ranked_features_tpr %>%
  left_join(feature_map, by = c("feature" = "symbol")) %>%
  dplyr::select(method, print_name, importance)

# Print results table
str_out <- ""
for(i in 1:nrow(temp)) {
  str_out <- paste0(str_out,
                    temp[i,]$method,
                    " & ",
                    temp[i,]$print_name,
                    " & ",
                    temp[i,]$importance,
                    " \\\\ \\hline \n")
}
writeLines(str_out, file.path("output", "scratch.txt"))

temp <- ranked_features_fpr %>%
  left_join(feature_map, by = c("feature" = "symbol")) %>%
  dplyr::select(method, print_name, importance)

# Print results table
str_out <- ""
for(i in 1:nrow(temp)) {
  str_out <- paste0(str_out,
                    temp[i,]$method,
                    " & ",
                    temp[i,]$print_name,
                    " & ",
                    temp[i,]$importance,
                    " \\\\ \\hline \n")
}
writeLines(str_out, file.path("output", "scratch.txt"))

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

ranked_features_tpr %>%
  group_by(feature) %>%
  summarize(mean(rho))

ranked_features_fpr %>%
  group_by(feature) %>%
  summarize(mean(rho))

