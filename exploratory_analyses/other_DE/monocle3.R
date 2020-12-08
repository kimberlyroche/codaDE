library(devtools)
library(reticulate) # Python in R
library(dplyr)

# Installation instructions here:
# https://cole-trapnell-lab.github.io/monocle3/docs/installation/

# Done
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                        'limma', 'S4Vectors', 'SingleCellExperiment',
#                        'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

# Done
# devtools::install_github('cole-trapnell-lab/leidenbase')

# This seems to try to update packages that underlie tidyverse and break them
# in the process (rlang, vctrs). Luckily it looks like these can be reinstalled,
# loaded explicitly, tidyverse loaded, and then you can retry this install.
# Done
# devtools::install_github('cole-trapnell-lab/monocle3')

library(monocle3)

# Ensure umap-learn 0.30 is installed
# Note: This ran into errors. I ended up installing outside R using pip.
# (1) Install UMAP via
#       pip install umap-learn --user
# (2) Install this specific version of numpy which fixes a weird permissions error
#     in Windows (https://stackoverflow.com/questions/64729944/runtimeerror-the-current-numpy-installation-fails-to-pass-a-sanity-check-due-to)
#       pip install numpy==1.19.3
# (3) Install louvain via
#       pip install louvain --user
# (4) Test via
#       python()
#       import umap
#       import louvain
# reticulate::py_install('umap-learn', pip = T, pip_ignore_installed = T)

import("louvain")

source("exploratory_analyses/other_DE/run_first.R")

# Note: The gene and sample annotation data.frames must have rownames or you'll get stupid errors
# during the clustering. Reference: https://github.com/cole-trapnell-lab/monocle3/issues/233
sample_sheet_df <- data.frame(condition = labels)
rownames(sample_sheet_df) <- paste0("cell", 1:(n_samples_condition*2))
head(sample_sheet_df)
feature_data_df <- data.frame(gene_short_name = paste0("gene", 1:n_genes))
rownames(feature_data_df) <- paste0("gene", 1:n_genes)
head(feature_data_df)

# Create a new CellDataSet object
# This code changes from Monocle -> Monocle2 -> Monocle3 so watch out!
data <- new_cell_data_set(count_table,
                         cell_metadata = sample_sheet_df,
                         gene_metadata = feature_data_df)

# Following steps from: https://cole-trapnell-lab.github.io/monocle3/docs/starting/

## Step 1: Normalize and pre-process the data
# How to choose num_dim?
# Safe to ignore warnings about too-large number of features selected?
# https://github.com/satijalab/seurat/issues/1249
data <- preprocess_cds(data, num_dim = 100)

## Step 2: Batch effect correction (skipped)

## Step 3: Reduce the dimensions using UMAP
data <- reduce_dimension(data)

## Step 4: Cluster the cells
data <- cluster_cells(data)
plot_cells(data) # optional argument color_cells_by = "condition" (etc.)

# With regression:
gene_fits <- fit_models(data, model_formula_str = "~ condition")
fit_coefs <- coefficient_table(gene_fits)
condition_terms <- fit_coefs %>% filter(term == "condition")
condition_terms <- condition_terms %>% mutate(q_value = p.adjust(p_value))
sig_genes <- condition_terms %>% filter(q_value < 0.05) %>% pull(gene_short_name)
# length(sig_genes)/n_genes
quick_dirty_rate_calcs(sig_genes, feature_data_df[sim$da_genes,])

# With graph autocorrelation: (WORKS POORLY!)
pr_test_res <- graph_test(data, neighbor_graph = "knn", cores = 1)
sig_genes <- row.names(subset(pr_test_res, q_value < 0.05))
# length(sig_genes)/n_genes
quick_dirty_rate_calcs(sig_genes, feature_data_df[sim$da_genes,])

