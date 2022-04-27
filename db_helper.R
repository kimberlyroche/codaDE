source("path_fix.R")

library(RSQLite)
library(codaDE)
library(tidyverse)
library(SingleCellExperiment)
library(scran)
library(scater)

betas <- function(data, groups) {
  count_table <- t(data)
  n_genes <- nrow(count_table)
  cell_metadata <- data.frame(condition = groups)
  rownames(cell_metadata) <- colnames(count_table)
  rownames(count_table) <- paste0("gene", 1:n_genes)
  colnames(count_table) <- paste0("cell", 1:ncol(count_table))
  
  sce <- SingleCellExperiment(assays = list(counts = count_table),
                              colData = cell_metadata)
  clusters <- as.numeric(as.factor(groups))
  if(min(clusters) == 0) {
    clusters <- clusters + 1
  }
  sce <- suppressWarnings(computeSumFactors(sce, clusters = clusters))
  sce <- logNormCounts(sce)
  
  colLabels(sce) <- factor(clusters)
  
  sf <- sizeFactors(sce)
  norm_counts <- counts(sce)
  for(i in 1:ncol(norm_counts)) {
    norm_counts[,i] <- norm_counts[,i]/sf[i]
  }
  groups <- factor(groups)
  est <- numeric(nrow(norm_counts))
  for(i in 1:nrow(norm_counts)) {
    est[i] <- mean(norm_counts[i,groups == levels(groups)[2]])/mean(norm_counts[i,groups == levels(groups)[1]])
  }
  
  est
}

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

res <- dbGetQuery(conn, "SELECT * FROM results WHERE METHOD='scran'")
uuids <- unique(res$UUID)
upd <- 0
for(u in 3181:length(uuids)) {
  if(u %% 100 == 0) {
    cat(paste0("Iteration ", u, " / ", length(uuids), "\n"))
  }
  uuid <- uuids[u]
  # Pull data and recalculate betas
  data_obj <- readRDS(file.path("output", "datasets", paste0(uuid, ".rds")))
  groups <- data_obj$simulation$groups
  data <- data_obj$simulation$observed_counts1
  est <- betas(data, groups)
  new_betas <- paste(est, collapse = ";")
  # Where these betas are in use, wipe out our "rates" estimates
  upd <- upd + dbExecute(conn, paste0("UPDATE results SET BETAS='", new_betas, "', TPR=NULL, FPR=NULL WHERE UUID='", uuid, "' AND METHOD='scran' AND PARTIAL_INFO=0 AND BASELINE_TYPE='oracle' AND OBSERVED_TYPE='relative_abundances' AND BETA=2"))
  # Otherwise just update for good measure
  upd <- upd + dbExecute(conn, paste0("UPDATE results SET BETAS='", new_betas, "' WHERE UUID='", uuid, "' AND METHOD='scran' AND PARTIAL_INFO=0 AND BASELINE_TYPE='oracle' AND OBSERVED_TYPE='relative_abundances' AND BETA!=2"))
}

cat(paste0("Updated ", upd, " rows\n"))

dbDisconnect(conn)

