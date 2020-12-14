library(Seurat)
library(DESeq2)

source("exploratory_analyses/other_DE/run_first.R")

rownames(count_table) <- paste0("gene", 1:n_genes)
colnames(count_table) <- paste0("cell", 1:(n_samples_condition*2))
head(count_table)

# Create Seurat object manually:
# https://learn.gencore.bio.nyu.edu/single-cell-rnaseq/loading-your-own-data-in-seurat-reanalyze-a-different-dataset/
data <- CreateSeuratObject(counts = count_table, project = "simulation", min.cells = 1, min.features = 10)
DefaultAssay(data) <- "RNA"

# Note: Some kind of implicit filtering seems to happen here. Need to figure that out.
for(method in c("wilcox", "DESeq2")) {
  markers <- FindMarkers(data, ident.1 = colnames(count_table)[labels == 0],
                         ident.2 = colnames(count_table)[labels == 1],
                         test.use = method,
                         assay = "RNA")
  # Adjusted p-values are Bonferroni corrected
  sig_genes <- rownames(markers)[markers$p_val_adj < 0.05]
  quick_dirty_rate_calcs(sig_genes, rownames(count_table)[sim$da_genes])
}

# Note: Calling MAST with a Seurat object seems to throw an error I can't resolve.
# Update: I think this is because MAST expects to find a logCounts object somewhere. That's got to be in an obvious slot.
