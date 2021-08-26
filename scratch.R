source("path_fix.R")

library(tidyverse)
library(codaDE)
library(RSQLite)
library(gridExtra)
library(RColorBrewer)
library(randomForest)

source("ggplot_fix.R")

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

p <- 1000
partial <- 0
reference <- "self"

cat(paste0("Evaluating P=", p, ", PARTIAL=", partial, ", REF=", reference, "\n"))
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
                               "WHERE P=", p, " ",
                               "AND PARTIAL_INFO=", partial, " ",
                               "AND BASELINE_TYPE='",reference,"' ",
                               "AND FC_ABSOLUTE <= 10 ",
                               "AND FC_ABSOLUTE >= 0.1;"))

dbDisconnect(conn)

# Strip "result-less" entries
res <- res %>%
  filter(!is.na(TPR) & !is.na(FPR))

res$FC_plot <- sapply(res$FC_ABSOLUTE, function(x) {
  if(x < 1) {
    1 / x
  } else {
    x
  }
})
res$FC_plot <- cut(res$FC_plot, breaks = c(1, 2, 5, Inf))
levels(res$FC_plot) <- c("low", "moderate", "high")

# Select as "bad" data sets, those with high FC, high TPR, and high FPR
bad_sets <- res %>%
  filter(METHOD == "DESeq2") %>%
  filter(FC_plot == "high") %>%
  filter(FPR > 0.3) %>%
  filter(TPR > 0.8) %>%
  select(-c("ORACLE_BASELINE", "CALLS"))

# Select as "good" data sets, those with high FC, high TPR, and low FPR
good_sets <- res %>%
  filter(METHOD == "DESeq2") %>%
  filter(FC_plot == "high") %>%
  filter(FPR < 0.05) %>%
  filter(TPR > 0.8) %>%
  select(-c("ORACLE_BASELINE", "CALLS"))

# Do we have a reasonable sample of each?
dim(bad_sets)
dim(good_sets)

# Is correlation any different between "bad" and "good" datasets?
# table(bad_sets$CORRP)
# table(good_sets$CORRP)

uuids <- c(bad_sets$UUID, good_sets$UUID)
utype <- c(rep("bad", nrow(bad_sets)), rep("good", nrow(good_sets)))

results <- NULL

for(i in 1:length(uuids)) {
  ref_obj <- readRDS(file.path("output", "datasets", paste0(uuids[i], ".rds")))$simulation
  ref_data <- ref_obj$abundances
  p <- ncol(ref_data)
  n <- floor(nrow(ref_data)/2)
  m1 <- mean(rowSums(ref_data)[1:n])
  m2 <- mean(rowSums(ref_data)[(n+1):(n*2)])
  fc <- m2 / m1
  # cat(paste0("FC: ", round(fc, 1), "\n"))
  
  majority_sign <- sign(m2 - m1)
  # cat(paste0("Majority sign: ", majority_sign, "\n"))
  
  # # Visualize
  # groups <- ref_obj$groups
  # subset_features <- sample(1:ncol(ref_data), size = 500)
  # mean_A <- colMeans(ref_data[groups == unique(groups)[1],])
  # mean_B <- colMeans(ref_data[groups == unique(groups)[2],])
  # plot_bipartite_graph(log(mean_A[subset_features] + 0.5),
  #                      log(mean_B[subset_features] + 0.5))

  # How many features change in the same direction as overall change?
  m1 <- colMeans(ref_data[groups == unique(groups)[1],])
  m2 <- colMeans(ref_data[groups == unique(groups)[2],])
  q1 <- sum(sign(m2 - m1) == majority_sign) / p
  # cat(paste0("Proportion of feature in same direction as overall: ", round(q1, 2), "\n"))
  
  # The top 10% of features account for what % of the majority change?
  change <- m2 - m1
  majority_change <- change[sign(change) == majority_sign]
  ranks <- order(majority_change, decreasing = TRUE)
  n <- length(ranks)
  n10 <- round(n/10)
  q2 <- sum(majority_change[ranks[1:n10]]) / sum(majority_change)
  # cat(paste0("Top 10% of features are doing this percent of 'work': ", round(q2, 2), "\n"))
  
  # What's the max differential for a single feature as a fraction of majority change?
  q3 <- max(majority_change) / sum(majority_change)
  # cat(paste0("Max differential: ", round(q3, 2), "\n"))
  
  # Proportion "stable" features (i.e. non-negligibly abundant with less than
  # 1.5-fold change)
  include_features <- m2 >= 1 & m1 >= 1
  feature_fc <- m2[include_features] / m1[include_features]
  feature_fc[feature_fc < 1] <- 1 / feature_fc[feature_fc < 1]
  q4 <- sum(feature_fc < 1.5) / length(feature_fc)
  
  results <- rbind(results,
                   data.frame(uuid = uuids[i],
                              utype = utype[i],
                              fc = fc,
                              majority_sign = majority_sign,
                              prop_overall = q1,
                              p_top_10 = q2,
                              max_diff = q3,
                              p_stable = q4))
}

cat("Saving results...\n")
saveRDS(results, "bad_good_results.rds")
