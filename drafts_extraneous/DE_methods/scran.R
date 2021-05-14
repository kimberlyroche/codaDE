library(scran)
library(scater)
library(edgeR)

source("exploratory_analyses/other_DE/run_first.R")

cell_metadata <- data.frame(condition = labels)

rownames(cell_metadata) <- colnames(count_table)
rownames(count_table) <- paste0("gene", 1:n_genes)
colnames(count_table) <- paste0("cell", 1:(n_samples_condition*2))
head(count_table)

sce <- SingleCellExperiment(assays = list(counts = count_table),
                            colData = cell_metadata)

# Note: We'll need to add QC stuff here ultimately.

# Methods below taken from this tutorial:
# https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html

clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = clusters) # Use a first-pass clustering
sce <- logNormCounts(sce)

# Model gene variance
# Note: This looks goofy as hell width two completely discrete clusters.
dec <- modelGeneVar(sce)
# plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
# curve(metadata(dec)$trend(x), col="blue", add=TRUE)

# Get top 10% of highly variable genes
top.hvgs <- getTopHVGs(dec, prop = 0.1)
# top.hvgs <- getTopHVGs(dec, fdr.threshold = 0.05)

# These HVGs inform PCA
sce <- runPCA(sce, subset_row = top.hvgs)

# Automatically choose some number of clusters
output <- getClusteredPCs(reducedDim(sce))
npcs <- metadata(output)$chosen
cat("Number of informative dimensions inferred:",npcs,"\n")
reducedDim(sce, "PCAsub") <- reducedDim(sce, "PCA")[,1:npcs,drop=FALSE]

# Build a graph based on this truncated dimensionality reduction
g <- buildSNNGraph(sce, use.dimred = "PCAsub")
cluster <- igraph::cluster_walktrap(g)$membership

# Assigning to the 'colLabels' of the 'sce'
colLabels(sce) <- factor(cluster)
table(colLabels(sce))

# Visualize the clustering
sce <- runTSNE(sce, dimred = "PCAsub")
plotTSNE(sce, colour_by = "label", text_by = "label")

# markers <- findMarkers(sce) # t-test by default
markers <- findMarkers(sce, test.type = "wilcox") # t-test by default
filter_vector <- markers[[1]][,1:3]$FDR < 0.05
sig_genes <- rownames(markers[[1]])[filter_vector]
quick_dirty_rate_calcs(sig_genes, rownames(count_table)[sim$da_genes])

