library(edgeR)

source("exploratory_analyses/other_DE/run_first.R")

rownames(cell_metadata) <- colnames(count_table)
rownames(count_table) <- paste0("gene", 1:n_genes)
colnames(count_table) <- paste0("cell", 1:(n_samples_condition*2))
head(count_table)

sce <- SingleCellExperiment(assays = list(counts = count_table),
                            colData = cell_metadata)

# Confirm the false positive effect of just using library size
sce <- suppressWarnings(computeSumFactors(sce)) # warnings are about assuming UMI counts (ok)
dge_obj <- convertTo(sce, type = "edgeR")
sig_genes <- quick_dirty_edgeR_call(dge_obj)
quick_dirty_rate_calcs(sig_genes, rownames(count_table)[sim$da_genes])

# Confirm TMM normalization sucks too
dge_obj <- DGEList(counts = count_table, group = factor(labels))
dge_obj <- calcNormFactors(dge_obj, method = "TMM")
sig_genes <- quick_dirty_edgeR_call(dge_obj)
quick_dirty_rate_calcs(sig_genes, rownames(count_table)[sim$da_genes])
