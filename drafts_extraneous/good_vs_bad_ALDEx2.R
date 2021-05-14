library(tidyverse)
library(ALDEx2)
library(driver)

# Explore good- and bad-performing cases in ALDEx2

make_NB_calls <- function(data, groups) {
  sapply(1:ncol(data), function(idx) { call_DA_NB(data, groups, idx) } )
}

make_ALDEx2_calls <- function(data, groups) {
  # Spike-in a count for any all-zero features
  count_table <- t(data)
  n_genes <- nrow(count_table)
  n <- ncol(count_table)/2
  spike_idx <- which(rowSums(count_table) == 0)
  for(idx in spike_idx) {
    count_table[idx,1:n] <- rmultinom(1, size = 1, prob = rep(1/n, n))[,1]
    count_table[idx,(n+1):(n*2)] <- rmultinom(1, size = 1, prob = rep(1/n, n))[,1]
  }
  rownames(count_table) <- 1:n_genes
  colnames(count_table) <- 1:(n*2)
  res <- aldex(count_table, groups) # alt: denom = "iqlr"
  # Interpreting the output:
  # we.ep - Expected P value of Welch's t test
  # we.eBH - Expected Benjamini-Hochberg corrected P value of Welch's t test
  # wi.ep - Expected P value of Wilcoxon rank test
  # wi.eBH - Expected Benjamini-Hochberg corrected P value of Wilcoxon test
  # kw.ep - Expected P value of Kruskal-Wallace test
  # kw.eBH - Expected Benjamini-Hochberg corrected P value of Kruskal-Wallace test
  calls <- res$we.eBH
  names(calls) <- rownames(res)
  clr.counts <- clr_array(data + 1, parts = 2)
  return(list(calls = calls, clr.counts = clr.counts))
}

make_calls <- function(sim_data) {
  oracle_calls <- make_NB_calls(sim_data$abundances, sim_data$groups)
  results <- make_ALDEx2_calls(sim_data$observed_counts1, sim_data$groups)
  return(list(calls = results$calls, oracle_calls = oracle_calls, clr.counts = results$clr.counts))
}

# These are saved simulations (of hundreds) that are very similar in terms of
# data characteristics but nonetheless have very different performance after
# scran normalization.
#
# These have similar performance in terms of the baseline model (NB GLM), which
# calls ~50% features differential in both cases.

output_dir <- file.path("output", "ALDEx2_debugging")
alpha <- 0.05

# ------------------------------------------------------------------------------
#   Explore FALSE POSITIVES in "bad" case
# ------------------------------------------------------------------------------

# sim_data <- readRDS(file.path(output_dir, "high_fpr.rds"))
sim_data <- readRDS(file.path(output_dir, "low_fpr.rds"))
n <- length(sim_data$groups)/2

# Make calls
all_calls <- make_calls(sim_data)
calls <- all_calls$calls
oracle_calls <- all_calls$oracle_calls

# Calculate rates
de <- calls < alpha
sim_de <- oracle_calls < alpha
TP <- sum(de & sim_de)
FP <- sum(de & !sim_de)
TN <- sum(!de & !sim_de)
FN <- sum(!de & sim_de)

cat("TPR:", round(TP / (TP + FN), 3), "\n")
cat("FPR:", round(FP / (FP + TN), 3), "\n")

# Proportion dropouts observed
dropout_incr <- sum(colSums(sim_data$observed_counts1[(n+1):(n*2),]) < 1) /
  sum(colSums(sim_data$observed_counts1[1:n,]) < 1)
cat("Fold increase in dropouts:", round(dropout_incr, 3), "\n")

# Plot random false positives
FP_idx <- which(de & !sim_de)
gene_idx <- sample(FP_idx, size = 1)

TN_idx <- which(!de & !sim_de)
gene_idx <- sample(TN_idx, size = 1)
plot_df <- data.frame(index = rep(1:(n*2), 3),
                      abundance = c(sim_data$abundances[,gene_idx],
                                    sim_data$observed_counts1[,gene_idx],
                                    all_calls$clr.counts[,gene_idx]),
                      type = c(rep("A", n*2),
                               rep("B", n*2),
                               rep("C", n*2)),
                      condition = rep(c(rep("A", n), rep("B", n)), 3))
plot_df$type <- factor(plot_df$type)
levels(plot_df$type) <- c("true abundances",
                          "resampled abundances",
                          "CLR resampled counts")
plot_df$condition <- factor(plot_df$condition)
ggplot(plot_df, aes(x = index, y = abundance, color = condition)) +
  geom_point(size = 3) +
  facet_wrap(. ~ type, scales = "free_y") +
  labs(x = "sample index",
       y = "(log) abundance")
ggsave(file.path(output_dir, "bad_case_FP.png"),
       units = "in",
       dpi = 100,
       height = 4,
       width = 12)
