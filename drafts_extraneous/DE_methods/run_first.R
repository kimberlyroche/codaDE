library(codaDE)
library(scran)
library(scater)

quick_dirty_edgeR_call <- function(dge_obj) {
  design <- model.matrix(~ colData(sce)$condition)
  dge_obj <- estimateDisp(dge_obj, design)
  # LRT recommended for single-cell data
  fit <- glmFit(dge_obj, design)
  lrt <- glmLRT(fit, coef = 2)
  is.de <- decideTestsDGE(lrt)
  rownames(is.de)[is.de != 0]
}

quick_dirty_rate_calcs <- function(test, reference) {
  tp <- sum(test %in% reference)
  fp <- sum(!(test %in% reference))
  fn <- sum(!(reference %in% test))
  tn <- n_genes - tp - fp - fn
  cat("TPR:", round(tp/(tp + fn), 2) ,"\n")
  cat("FPR:", round(fp/(fp + tn), 2) ,"\n")
}

# # Generate dummy count data -- copied from monocle3.R
n_genes <- 5000
n_samples_condition <- 100
# dispersion <- 0.1
# mean_multiplier <- 1.5
# count_table <- matrix(rnbinom(n_genes*n_samples_condition*2, mu = 100, size = 1/dispersion), n_genes, n_samples_condition*2)
# # Simulate very dumb DE for the first 1/2 of the genes in the second condition
# count_table[1:(n_genes/2),(n_samples_condition+1):(n_samples_condition*2)] <- rnbinom((n_genes/2)*n_samples_condition, mu = 100*mean_multiplier, size = 1/dispersion)
# # Poor man's TPM
# count_table <- apply(count_table, 2, function(x) round((x/sum(x))*1e6))
# cell_metadata <- data.frame(batch = c(rep("A", n_samples_condition), rep("B", n_samples_condition)))

# Better simulation
sim <- simulate_singlecell_RNAseq(p = n_genes, n = n_samples_condition, k = 1, proportion_da = 0.6, size_factor_correlation = 0)
true_counts <- t(sim$abundances)
count_table <- t(sim$observed_counts)
labels <- sim$groups
# cell_metadata <- data.frame(batch = labels)
# 
# rownames(cell_metadata) <- colnames(count_table)
# rownames(count_table) <- paste0("gene", 1:n_genes)
# colnames(count_table) <- paste0("cell", 1:(n_samples_condition*2))
# head(count_table)
# 
# sce <- SingleCellExperiment(assays = list(counts = count_table),
#                             colData = cell_metadata)

# Confirm the false positive effect of just using library size
# sce <- suppressWarnings(computeSumFactors(sce)) # warnings are about assuming UMI counts (ok)
# dge_obj <- convertTo(sce, type = "edgeR")
# quick_dirty_edgeR_call(dge_obj)

# Confirm TMM normalization sucks too
# dge_obj <- DGEList(counts = count_table, group = factor(labels))
# dge_obj <- calcNormFactors(dge_obj, method = "TMM")
# quick_dirty_edgeR_call(dge_obj)

# Confirm that cluster size factor normalization ("pooling and deconvolution") works well
# clusters <- quickCluster(sce)
# sce <- computeSumFactors(sce, clusters = clusters)
# dge_obj <- convertTo(sce, type = "edgeR")
# quick_dirty_edgeR_call(dge_obj)







