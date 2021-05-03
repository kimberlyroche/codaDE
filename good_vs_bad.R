library(codaDE)
library(tidyverse)
library(gridArrange)
library(uuid)
library(SingleCellExperiment)
library(scran)
library(scater)

p <- 100
correlation_prop <- 0

base_correlation <- diag(p)
concentration <- 1e6

n <- 10 # replicate number
asymmetry <- 1
proportion_da <- 0.75
spike_in <- FALSE
iterations <- 400

palette <- generate_highcontrast_palette(p)

# Copied from iterate_jobs.R

DE_by_NB <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- sapply(1:p, function(idx) { call_DA_NB(ref_data, groups, idx) } )
  }
  calls <- sapply(1:p, function(idx) { call_DA_NB(data, groups, idx) } )
  return(list(oracle_calls = oracle_calls, calls = calls))
}

DE_by_scran <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_scran(ref_data, groups)
  }
  calls <- call_DA_scran(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

DE_by_ALDEx2 <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_ALDEx2(ref_data, groups)
  }
  calls <- call_DA_ALDEx2(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

calc_DE_discrepancy <- function(ref_data, data, groups, method = "NB",
                                oracle_calls = NULL) {
  if(method == "scran") {
    DE_calls <- DE_by_scran(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls
    }
    calls <- DE_calls$calls
  } else if(method == "ALDEx2") {
    DE_calls <- DE_by_ALDEx2(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls
    }
    calls <- DE_calls$calls
  } else {
    DE_calls <- DE_by_NB(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls
    }
    calls <- DE_calls$calls
  }
  
  de <- calls < 0.05
  sim_de <- oracle_calls < 0.05
  
  TP <- sum(de & sim_de)
  FP <- sum(de & !sim_de)
  
  TN <- sum(!de & !sim_de)
  FN <- sum(!de & sim_de)
  
  return(list(tpr = TP/(TP+FN), fpr = FP/(FP+TN), oracle_calls = oracle_calls))  
}

k_percent <- function(abundances, percent = 0.9) {
  cond1 <- abundances[1:n,]
  cond2 <- abundances[(n+1):(n*2),]
  delta_features <- abs(colMeans(cond2) - colMeans(cond1))
  rank_da <- order(delta_features, decreasing = TRUE)
  rank_deltas <- delta_features[rank_da]
  rank_deltas <- rank_deltas / sum(rank_deltas)
  k_prop <- rank_deltas[1]
  for(i in 2:p) {
    if(k_prop > percent) {
      break
    }
    k_prop <- k_prop + rank_deltas[i]
  }
  return(c(i, k_prop))
}

data_file <- file.path("output", "temp", "data.rds")
if(file.exists(data_file)) {
  data <- readRDS(data_file)
} else {
  data <- data.frame(lfc = c(),
                     tpr = c(),
                     fpr = c(),
                     uuid = c())
  for(i in 1:iterations) {
    data_obj <- build_simulated_reference(p = p,
                                          log_mean = 1,
                                          log_var = 2,
                                          log_noise_var = 2,
                                          base_correlation = base_correlation,
                                          concentration = concentration)
    sim_data <- simulate_sequence_counts(n = n,
                                          p = p,
                                          data_obj = data_obj,
                                          asymmetry = asymmetry,
                                          proportion_da = proportion_da,
                                          spike_in = spike_in)
    r1 <- rowSums(sim_data$abundances[1:n,])
    r2 <- rowSums(sim_data$abundances[(n+1):(n*2),])
    m1 <- mean(r1)
    m2 <- mean(r2)
    delta_mean_v2 <- max(c(m1, m2)) / min(c(m1, m2))
    rates_baseline <- calc_DE_discrepancy(sim_data$abundances[,1:p], sim_data$observed_counts1[,1:p], sim_data$groups)
    oracle_calls <- rates_baseline$oracle_calls
    rates <- calc_DE_discrepancy(sim_data$abundances[,1:p],
                                 sim_data$observed_counts1[,1:p],
                                 sim_data$groups,
                                 method = "scran",
                                 oracle_calls = oracle_calls)
    uuid <- UUIDgenerate()
    data <- rbind(data,
                  data.frame(lfc = log(delta_mean_v2),
                       tpr = rates$tpr,
                       fpr = rates$fpr,
                       uuid = uuid))
    cat(nrow(data), "\n")
    saveRDS(sim_data, file = file.path("output", "temp", paste0(uuid, ".rds")))
  }
  data <- cbind(data, id = 1:nrow(data))
  saveRDS(data, file = data_file)
}

pl1 <- ggplot(data, aes(x = fpr, y = tpr, color = lfc)) +
  geom_point(size = 4) +
  scale_color_gradient2(low = "blue", mid = "#BBBBBB", high = "red", midpoint = mean(data$lfc))
pl2 <- ggplot(data, aes(x = fpr, y = tpr, label = id)) +
  geom_text(color = "black", size = 3)
grid.arrange(grobs = list(pl1, pl2), ncol = 2)

# ------------------------------------------------------------------------------
#   Color plot for number of features accounting for ~90% of DE
#   (This is not an obvious driver of differences in FPR)
# ------------------------------------------------------------------------------

k <- numeric(nrow(data))
for(i in 1:iterations) {
  sim_data <- readRDS(file.path("output", "temp", paste0(data$uuid[i], ".rds")))
  k[i] <- k_percent(sim_data$abundances)[1]
  # r1 <- rowSums(sim_data$abundances[1:n,])
  # r2 <- rowSums(sim_data$abundances[(n+1):(n*2),])
  # m1 <- mean(r1)
  # m2 <- mean(r2)
  # k[i] <- log(m2 / m1)
}

plot_data <- cbind(data, k = k)
ggplot(plot_data, aes(x = fpr, y = tpr, color = k)) +
  geom_point(size = 4) +
  scale_color_gradient2(low = "blue",
                        mid = "#BBBBBB",
                        high = "red",
                        midpoint = mean(plot_data$k))

# ------------------------------------------------------------------------------
#   Evaluate individual simulations
# ------------------------------------------------------------------------------

idx <- sample(which(data$fpr < 0.05), size = 1) # good
# idx <- sample(which(data$fpr > 0.15), size = 1) # bad

sim_data <- readRDS(file.path("output", "temp", paste0(data[data$id == idx,]$uuid, ".rds")))

# count_table <- as.data.frame(t(sim_data$abundances)) # features x samples
# count_table <- cbind(1:nrow(count_table), count_table)
# colnames(count_table) <- c("feature", 1:(ncol(count_table)-1))
# counts_long <- pivot_longer(count_table, !feature, names_to = "sample", values_to = "abundance")
# counts_long$feature <- factor(counts_long$feature)
# counts_long$sample <- factor(counts_long$sample, levels = 1:(n*2))
# plot_stacked_bars(counts_long, palette = palette, save_name = NULL)

# rates_baseline <- calc_DE_discrepancy(sim_data$abundances[,1:p], sim_data$observed_counts1[,1:p], sim_data$groups)
# oracle_calls <- rates_baseline$oracle_calls
# DE_calls <- DE_by_scran(sim_data$abundances[,1:p], sim_data$observed_counts1[,1:p], sim_data$groups, oracle_calls = oracle_calls)

count_table <- t(sim_data$observed_counts1[,1:p])
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
# sce <- runTSNE(sce, dimred = "PCAsub")
# plotTSNE(sce, colour_by = "label", text_by = "label")
plotPCA(sce, colour_by = "label", text_by = "label")


