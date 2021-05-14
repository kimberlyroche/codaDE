library(tidyverse)
library(SingleCellExperiment)
library(scran)
library(scater)

# Explore good- and bad-performing cases in scran.

make_NB_calls <- function(data, groups) {
  sapply(1:ncol(data), function(idx) { call_DA_NB(data, groups, idx) } )
}

make_scran_calls <- function(data, groups) {
  count_table <- t(data)
  n_genes <- nrow(count_table)
  n_samples_condition <- ncol(count_table)/2
  cell_metadata <- data.frame(condition = groups)
  rownames(cell_metadata) <- colnames(count_table)
  rownames(count_table) <- paste0("gene", 1:n_genes)
  colnames(count_table) <- paste0("cell", 1:(n_samples_condition*2))
  sce <- SingleCellExperiment(assays = list(counts = count_table),
                              colData = cell_metadata)
  clusters <- suppressWarnings(quickCluster(sce, min.size = 1))
  sce <- suppressWarnings(computeSumFactors(sce, clusters = clusters)) # Use a first-pass clustering
  sce <- logNormCounts(sce)
  sce <- suppressWarnings(scater::runPCA(sce))
  output <- getClusteredPCs(reducedDim(sce))
  npcs <- metadata(output)$chosen # number of "informative dimensions"
  reducedDim(sce, "PCAsub") <- reducedDim(sce, "PCA")[,1:npcs,drop=FALSE]
  g <- buildSNNGraph(sce, use.dimred = "PCAsub")
  cluster <- igraph::cluster_walktrap(g)$membership
  colLabels(sce) <- factor(cluster)
  log_expr <- assay(sce, "logcounts")
  
  markers <- findMarkers(sce, test.type = "wilcox") # t-test by default
  filter_vector <- markers[[1]][,1:3]$FDR < 0.05
  sig_genes <- rownames(markers[[1]])[filter_vector]
  calls <- rep(1, n_genes)
  calls[which(rownames(count_table) %in% sig_genes)] <- 0
  return(list(calls = calls, log_expr = log_expr))
}

make_calls <- function(sim_data) {
  oracle_calls <- make_NB_calls(sim_data$abundances, sim_data$groups)
  results <- make_scran_calls(sim_data$observed_counts1, sim_data$groups)
  return(list(calls = results$calls, oracle_calls = oracle_calls, scran_log_expr = results$log_expr))
}

# These are saved simulations (of hundreds) that are very similar in terms of
# data characteristics but nonetheless have very different performance after
# scran normalization.
#
# These have similar performance in terms of the baseline model (NB GLM), which
# calls ~50% features differential in both cases.

output_dir <- file.path("output", "scran_debugging")
alpha <- 0.05

# ------------------------------------------------------------------------------
#   Explore FALSE POSITIVES in "bad" case
# ------------------------------------------------------------------------------

sim_data <- readRDS(file.path(output_dir, "high_fpr.rds"))
# sim_data <- readRDS(file.path(output_dir, "low_fpr.rds"))
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
plot_df <- data.frame(index = rep(1:(n*2), 3),
                      abundance = c(sim_data$abundances[,gene_idx],
                                    sim_data$observed_counts1[,gene_idx],
                                    all_calls$scran_log_expr[gene_idx,]),
                      type = c(rep("A", n*2),
                               rep("B", n*2),
                               rep("C", n*2)),
                      condition = rep(c(rep("A", n), rep("B", n)), 3))
plot_df$type <- factor(plot_df$type)
levels(plot_df$type) <- c("true abundances",
                          "resampled abundances",
                          "scran corrected log counts")
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
