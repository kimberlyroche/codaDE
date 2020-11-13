library(Matrix)
library(ggplot2)
library(stringr)
library(MASS)
library(edgeR)

## --------------------------------------------------------------------------------------------------------------------------------
##   Some notes
##     - Conversion between data representations (e.g. TPM) described in this paper:
##         https://rnajournal.cshlp.org/content/early/2020/04/13/rna.074922.120.full.pdf
##     - I'm visually validating cell types by looking at levels of insulin (Ins2), glucagon (Gcg), and somatostatin (Sst). This
##       seems like standard pracice. Per DiGruccio et al. paper:
##         "Cells that had an RPKM value > 10 k of either Sst, or Ins2, or Gcg were defined as delta, beta, or alpha cells,
##         respectively."
## --------------------------------------------------------------------------------------------------------------------------------

#setwd("C:/Users/kim/Documents/codaDE/pancreatic_datasets")
setwd("/data/mukherjeelab/roche/codaDE/exploratory/pancreatic_datasets")

bulk_dir <- "Adriaenssens2016"
sc_dir <- "MarquinaSanchez2020"

# These are quality control settings used by Marquina-Sanchez et al.
MIN.UMIS <- 500
MIN.GENES <- 200
SCALE.FACTOR <- 10^6

show_images <- TRUE

## --------------------------------------------------------------------------------------------------------------------------------
##   FUNCTIONS
## --------------------------------------------------------------------------------------------------------------------------------

# Assumes gene_idx is 0-indexed and sample_IDs is a set of cell barcodes.
# pull_expr_vector <- function(mat_obj, gene_idx, sample_IDs = NULL) {
#   sample_idx <- 1:ncol(mat_obj)
#   if(!is.null(sample_IDs)) {
#     sample_idx <- sort(which(colnames(sc) %in% sample_IDs))
#   }
#   # Remember to switch indexing back to 1-indexing for this accessor!
#   mat_obj[gene_idx + 1, sample_idx]
# }

umi_to_tpm <- function(mat_obj) {
  tpm.sc <- t(t(mat_obj) * SCALE.FACTOR/Matrix::colSums(mat_obj))
  as(tpm.sc, "dgTMatrix")
}

pull_bulk_expr <- function(bulk_expr, gene_idx) {
  bulk_alpha <- round(unname(bulk_expr[gene_idx, colnames(bulk_expr) == "Alpha"]))
  bulk_beta <- round(unname(bulk_expr[gene_idx, colnames(bulk_expr) == "Beta"]))
  list(alpha = bulk_alpha, beta = bulk_beta)
}

pull_sc_expr <- function(sc_expr, gene_idx, sc_alpha_samples, sc_beta_samples) {
  sc_alpha <- round(sc_expr[gene_idx, sc_alpha_samples])
  sc_beta <- round(sc_expr[gene_idx, sc_beta_samples])
  list(alpha = sc_alpha, beta = sc_beta)
}

quick_plot_expr <- function(bulk_expr_obj, sc_expr_obj, named_gene, sc_alpha_idx, sc_beta_idx) {
  gene_name <- names(named_gene)
  gene_idx <- unname(named_gene)
  bulk_expr_vector <- pull_bulk_expr(bulk_expr_obj, gene_idx)
  bulk_celltypes <- c(rep("alpha", length(bulk_expr_vector$alpha)), rep("beta", length(bulk_expr_vector$beta)))
  sc_expr_vector <- pull_sc_expr(sc_expr_obj, gene_idx, sc_alpha_idx, sc_beta_idx)
  sc_celltypes <- character(ncol(sc_expr_obj))
  sc_celltypes[sc_alpha_idx] <- "alpha"
  sc_celltypes[sc_beta_idx] <- "beta"
  df <- data.frame(index = 1:length(bulk_celltypes),
                   expression = c(bulk_expr_vector$alpha, bulk_expr_vector$beta),
                   type = "bulk",
                   celltype = factor(bulk_celltypes))
  df <- rbind(df, data.frame(index = 1:length(sc_celltypes),
                             expression = c(sc_expr_vector$alpha, sc_expr_vector$beta),
                             type = "SC",
                             celltype = factor(sc_celltypes)))
  ggplot(df, aes(x = index, y = expression, color = celltype)) +
    geom_point() +
    facet_wrap(~ type, ncol = 2, scales = "free")
}

# call_set is TP, FP, TN, or FN
visualize_calls <- function(bulk_expr_obj, sc_expr_obj, call_set, sc_alpha_idx, sc_beta_idx) {
  gene <- sample(call_set)[1]
  gene_obj <- which(rownames(bulk_subset) == gene)
  names(gene_obj) <- gene
  quick_plot_expr(bulk_subset, sc_subset, gene_obj, 1:100, 101:200)
}

## --------------------------------------------------------------------------------------------------------------------------------
##   PARSE BULK DATA
## --------------------------------------------------------------------------------------------------------------------------------

bulk <- read.table(file.path(bulk_dir, "GSE76017_data_geo.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# The bulk data has "duplicate" gene names. We'll average these for now but we should probably take a closer look at this.
# I assume they're isoforms of the same gene. Maybe we want to think of another way of dealing with these.
# The below takes ~2 min. to run.
unique_genes <- sort(unique(bulk$gene_short_name))
bulk.collapsed <- matrix(NA, length(unique_genes), ncol(bulk)-2)
for(ug in 1:length(unique_genes)) {
  gene_name <- unique_genes[ug]
  selected_rows <- bulk[bulk$gene_short_name == gene_name,3:ncol(bulk)]
  collapsed_rows <- colMeans(selected_rows)
  bulk.collapsed[ug,] <- collapsed_rows
}
rownames(bulk.collapsed) <- unique_genes
bulk_celltypes <- unname(sapply(colnames(bulk)[3:ncol(bulk)], function(x) {
  str_match(x, "^([[:alpha:]]+)\\.")[1,2]
}))
colnames(bulk.collapsed) <- bulk_celltypes
bulk <- bulk.collapsed
rm(bulk.collapsed)

sst_idx_bulk <- which(rownames(bulk) == "Sst") # delta
ins2_idx_bulk <- which(rownames(bulk) == "Ins2") # beta
gcg_idx_bulk <- which(rownames(bulk) == "Gcg") # alpha
gapdh_idx_bulk <- which(rownames(bulk) == "Gapdh") # housekeeping

df <- data.frame(expression = unlist(bulk[sst_idx_bulk,]), gene = "Sst", celltype = bulk_celltypes)
df <- rbind(df, data.frame(expression = unlist(bulk[ins2_idx_bulk,]), gene = "Ins2", celltype = bulk_celltypes))
df <- rbind(df, data.frame(expression = unlist(bulk[gcg_idx_bulk,]), gene = "Gcg", celltype = bulk_celltypes))
df <- rbind(df, data.frame(expression = unlist(bulk[gapdh_idx_bulk,]), gene = "Gapdh", celltype = bulk_celltypes))
p <- ggplot(df, aes(x = as.factor(celltype), y = expression)) +
  geom_boxplot() +
  facet_wrap(~ gene, scales = "free") +
  xlab("Cell type") +
  ylab("Bulk reads")
if(show_images) {
  show(p)
} else {
  ggsave("01_marker_expression_in_bulk.png", p, units = "in", dpi = 100, height = 4, width = 6)
}

## --------------------------------------------------------------------------------------------------------------------------------
##   PARSE SINGLE-CELL DATA
## --------------------------------------------------------------------------------------------------------------------------------

# Parse the annotation file.
sc_annot <- read.table(file.path(sc_dir, "GSE142465_MouseLTI_CellAnnotation_final.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# Some stats:
#   - Annotations include `celltype` and `celltype2`. These agree in 92% of cells. I think they're generated by two sets of
#     heuristics. I'm arbitrarily using cells in the intersection of these two labels.
#   - Approx. 97% of sample IDs (barcodes) are shared between the count data and annotation files. The paper seems to imply the
#     rest might have been contaminated, very poor quality, etc.

# The below code was copied from 10X's website:
#   https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices
# This parses the pieces into a dgTMatrix sparse matrix object (from package `Matrix`) and should be basically equivalent to the
# Read10X function in Seurat that goes: sc <- Read10X(data.dir = "MarquinaSanchez2020/10X", gene.column = 1)
# This code takes ~2 min. to run.
barcode.path <- file.path(sc_dir, "GSE142465_MouseLTI_correctedTPM_Matrix_barcodes.txt.gz")
features.path <- file.path(sc_dir, "GSE142465_MouseLTI_correctedTPM_Matrix_genes.txt.gz")
matrix.path <- file.path(sc_dir, "GSE142465_MouseLTI_correctedTPM_Matrix.mtx.gz")
sc <- readMM(file = matrix.path)
feature.names <- read.delim(features.path, 
                            header = FALSE,
                            stringsAsFactors = FALSE)
barcode.names <- read.delim(barcode.path, 
                            header = FALSE,
                            stringsAsFactors = FALSE)
colnames(sc) <- barcode.names$V1 # sample IDs
rownames(sc) <- feature.names$V1 # gene IDs
# temp <- sc

# This dgTMatrix matrix works like so: (i,j) indexes the row of the ith gene and the jth sample and x is
# the value of that index. Missing (i,j) indices are implicitly zero!
# The alternative dgCMatrix matrix format uses (i,p) indexing: i indexes the row of the ith gene and p is
# an offset, capturing the number of non-zero entries in a given column, e.g. p[1] = 0 always, p[2] = number
# of non-zero entries in column 1, etc.
# Note that both use 0-indexing and can be translated between each other with `as(obj, "dgCMatrix")` etc.

# Filter down to the control (DMSO) samples via metadata.
sc_annot <- sc_annot[which(sc_annot$treatment == "DMSO"),]

# Filter down to consensus alpha, beta, and delta cells via metadata.
sc_annot <- sc_annot[which(sc_annot$celltype == sc_annot$celltype2),]
sc_annot <- sc_annot[which(sc_annot$celltype %in% c("Alpha", "Beta")),]

# Quality filter as in the paper source code. (See 05_PrepSpikeIns_Clean.R.) Everything identified above (DMSO x cell type) seems
# to pass QC.
sc <- sc[,colnames(sc) %in% sc_annot$rn]
cell.umis <- Matrix::colSums(sc)
sc <- sc[,cell.umis >= MIN.UMIS]
cell.genes <- Matrix::colSums(sc != 0)
sc <- sc[,cell.genes >= MIN.GENES]

# Get the sample IDs in the annotations ordered as in the count data file for easier indexing.
sc_annot <- sc_annot[order(sc_annot$rn),]
sc <- sc[,order(colnames(sc))]

sample_celltypes <- sc_annot$celltype
n_samples_sc <- ncol(sc)
n_genes_sc <- nrow(sc)

ins2_idx_sc <- which(rownames(sc) == "Ins2") # 1-indexed
gcg_idx_sc <- which(rownames(sc) == "Gcg")
gapdh_idx_sc <- which(rownames(sc) == "Gapdh")
eef2_idx_sc <- which(rownames(sc) == "Eef2")

df <- data.frame(expression = sc[ins2_idx_sc,], gene = "Ins2", celltype = sample_celltypes)
df <- rbind(df, data.frame(expression = sc[gcg_idx_sc,], gene = "Gcg", celltype = sample_celltypes))
df <- rbind(df, data.frame(expression = sc[gapdh_idx_sc,], gene = "Gapdh", celltype = sample_celltypes))
df <- rbind(df, data.frame(expression = sc[eef2_idx_sc,], gene = "Eef2", celltype = sample_celltypes))
p <- ggplot(df, aes(x = as.factor(celltype), y = expression)) +
  geom_jitter(width = 0.2) +
  facet_wrap(~ gene, scales = "free") +
  xlab("Cell type") +
  ylab("UMI counts")
if(show_images) {
  show(p)
} else {
  ggsave("02_marker_expression_in_sc.png", p, units = "in", dpi = 100, height = 4, width = 6)
}

# The total UMIs is systematically lower in alpha cells.
alpha_idx <- which(sc_annot$celltype == "Alpha") # 1-indexed
beta_idx <- which(sc_annot$celltype == "Beta")

df <- data.frame(nUMI = Matrix::colSums(sc[,alpha_idx]), celltype = "Alpha")
df <- rbind(df, data.frame(nUMI = Matrix::colSums(sc[,beta_idx]), celltype = "Beta"))
p <- ggplot(df, aes(x = as.factor(celltype), y = nUMI)) +
  geom_jitter(width = 0.2) +
  xlab("Cell type") +
  ylab("UMI counts")
if(show_images) {
  show(p)
} else {
  ggsave("03_marker_expression_in_sc.png", p, units = "in", dpi = 100, height = 4, width = 6)
}

# Repeat the first visualization with TPM-normalized data.
sc_tpm <- umi_to_tpm(sc)
df <- data.frame(expression = sc_tpm[ins2_idx_sc,], gene = "Ins2", celltype = sample_celltypes)
df <- rbind(df, data.frame(expression = sc_tpm[gcg_idx_sc,], gene = "Gcg", celltype = sample_celltypes))
df <- rbind(df, data.frame(expression = sc_tpm[gapdh_idx_sc,], gene = "Gapdh", celltype = sample_celltypes))
df <- rbind(df, data.frame(expression = sc_tpm[eef2_idx_sc,], gene = "Eef2", celltype = sample_celltypes))
p <- ggplot(df, aes(x = as.factor(celltype), y = expression)) +
  geom_jitter(width = 0.2) +
  facet_wrap(~ gene, scales = "free") +
  xlab("Cell type") +
  ylab("UMI counts")
if(show_images) {
  show(p)
} else {
  ggsave("04_marker_expression_in_sc_TPM.png", p, units = "in", dpi = 100, height = 4, width = 6)
}

# Gapdh x library size (nUMI)
df <- data.frame(x = Matrix::colSums(sc), y = sc[gapdh_idx_sc,])
p <- ggplot(df, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  xlab("nUMI") +
  ylab("Gadph expression")
if(show_images) {
  show(p)
} else {
  ggsave("05_marker_expression_in_sc_TPM.png", p, units = "in", dpi = 100, height = 6, width = 6)
}

## --------------------------------------------------------------------------------------------------------------------------------
##   PLOTTING BULK VS. SINGLE-CELL
## --------------------------------------------------------------------------------------------------------------------------------

# Get intersection of genes between bulk and single-cell and chop down each data set to this intersection. We'll need to
# assess what the missing genes are becuase I'm making the assumption the gene spaces are the same here (!).
shared_genes <- sort(intersect(rownames(bulk), rownames(sc)))
n_shared <- length(shared_genes)
all_genes <- union(rownames(bulk), rownames(sc))
cat("Proportion of genes in intersection vs. union:",round(n_shared/length(all_genes), 2),"\n")
cat("Number of genes in common:",n_shared,"\n")

# Order genes (rows) equivalently in intersected bulk and single-cell data.
sc_shared <- sc[rownames(sc) %in% shared_genes,]
sc_shared <- sc_shared[order(rownames(sc_shared)),]
bulk_shared <- bulk[rownames(bulk) %in% shared_genes,]
bulk_shared <- bulk_shared[order(rownames(bulk_shared)),]

# Get mean expression for genes in each data set.
sc_expr_means <- Matrix::rowMeans(sc_shared) # nUMI
bulk_expr_means <- rowMeans(bulk_shared)

quantile_each <- 0.95
plot_idx <- sc_expr_means <= quantile(sc_expr_means, probs = c(quantile_each)) & bulk_expr_means <= quantile(bulk_expr_means, probs = c(quantile_each))
df <- data.frame(x = bulk_expr_means[plot_idx], y = sc_expr_means[plot_idx])
p <- ggplot(df, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  xlab("Mean bulk expression (counts)") +
  ylab("Mean single-cell expression (UMI)")
if(show_images) {
  show(p)
} else {
  ggsave("06_bulk_vs_sc.png", p, units = "in", dpi = 100, height = 6, width = 6)
}

# Now do the TPM x TPM version.
sc_tpm <- umi_to_tpm(sc_shared)
bulk_tpm <- t(t(bulk_shared) * SCALE.FACTOR/colSums(bulk_shared))

sc_expr_means <- Matrix::rowMeans(sc_tpm)
bulk_expr_means <- rowMeans(bulk_tpm)

quantile_each <- 0.95
plot_idx <- sc_expr_means <= quantile(sc_expr_means, probs = c(quantile_each)) & bulk_expr_means <= quantile(bulk_expr_means, probs = c(quantile_each))
df <- data.frame(x = bulk_expr_means[plot_idx], y = sc_expr_means[plot_idx])
p <- ggplot(df, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  xlab("Mean bulk expression (counts)") +
  ylab("Mean single-cell expression (UMI)")
if(show_images) {
  show(p)
} else {
  ggsave("07_bulk_vs_sc_TPM.png", p, units = "in", dpi = 100, height = 6, width = 6)
}

## --------------------------------------------------------------------------------------------------------------------------------
##   DIFFERENTIAL EXPRESSION CALLING WITH edgeR
## --------------------------------------------------------------------------------------------------------------------------------

alpha <- 0.05

sc_beta_samples <- which(sc_annot$celltype == "Beta")
sc_alpha_samples <- which(sc_annot$celltype == "Alpha")

n_genes_to_test <- nrow(bulk_shared)
bulk_subset <- bulk_shared[1:n_genes_to_test,]
group_bulk <- colnames(bulk_subset)
group_sc <- c(sc_alpha_samples[1:100], sc_beta_samples[1:100])
sc_subset <- sc_shared[1:n_genes_to_test,group_sc]

de_calls_bulk <- numeric(n_genes_to_test)
de_calls_sc <- numeric(n_genes_to_test)

# From the edgeR quick start guide here:
# https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# Note: Column names must be unique for DGEList().

colnames(bulk_subset) <- paste0("sample", 1:ncol(bulk_subset))
y_bulk <- DGEList(counts = bulk_subset, group = factor(group_bulk))
keep_bulk <- filterByExpr(y_bulk)

colnames(sc_subset) <- paste0("sample", 1:ncol(sc_subset))
y_sc <- DGEList(counts = sc_subset, group = factor(group_sc))
keep_sc <- filterByExpr(y_sc)

keep_intersection <- keep_bulk & keep_sc
kept_genes <- names(keep_intersection)[keep_intersection]

y_bulk <- calcNormFactors(y_bulk)
design <- model.matrix(~ group_bulk)
y_bulk <- estimateDisp(y_bulk, design)
# quasi-likelihood recommended for bulk RNA-seq
fit <- glmQLFit(y_bulk, design)
qlf <- glmQLFTest(fit, coef = 2)
hits <- topTags(qlf, n = n_genes_to_test, p.value = 0.05)
hits_genes_bulk <- rownames(hits@.Data[[1]])
stable_genes_bulk <- kept_genes[!(kept_genes %in% hits_genes_bulk)]

y_sc <- calcNormFactors(y_sc)
design <- model.matrix(~ group_sc)
y_sc <- estimateDisp(y_sc, design)
# quasi-likelihood recommended for bulk RNA-seq
fit <- glmQLFit(y_sc, design)
qlf <- glmQLFTest(fit, coef = 2)
hits <- topTags(qlf, n = n_genes_to_test, p.value = 0.05)
hits_genes_sc <- rownames(hits@.Data[[1]])
stable_genes_sc <- kept_genes[!(kept_genes %in% hits_genes_sc)]

# Bulk RNA-seq calls are "truth"
TP <- intersect(hits_genes_bulk, hits_genes_sc)
FN <- intersect(hits_genes_bulk, stable_genes_sc)
TN <- intersect(stable_genes_bulk, stable_genes_sc)
FP <- intersect(stable_genes_bulk, hits_genes_sc)

cat("Genes evaluated:",length(kept_genes),"/",n_genes_to_test,"\n")
cat("TP:",length(TP),"/",length(kept_genes),"/",round(length(TP)/length(kept_genes), 2)*100,"\n")
cat("FN:",length(FN),"/",length(kept_genes),"/",round(length(FN)/length(kept_genes), 2)*100,"\n")
cat("TN:",length(TN),"/",length(kept_genes),"/",round(length(TN)/length(kept_genes), 2)*100,"\n")
cat("FP:",length(FP),"/",length(kept_genes),"/",round(length(FP)/length(kept_genes), 2)*100,"\n")

# Restore column labels to bulk_subset for subsequent visualization.
colnames(bulk_subset) <- group_bulk

# Which genes were excluded by "keep"? Ans: These seem to be genes with constant near-zero expression.
par(mfrow = c(1,2))

kept <- sample(which(keep_bulk == TRUE))[1]
quick_plot_expr(bulk_subset, kept)
quick_plot_expr(sc_subset, kept, sc_alpha_samples = 1:100, sc_beta_samples = 101:200)

discarded <- sample(which(keep_bulk == FALSE))[1]
quick_plot_expr(bulk_subset, discarded)
quick_plot_expr(sc_subset, discarded, sc_alpha_samples = 1:100, sc_beta_samples = 101:200)

# What do hits and misses look like?

p <- visualize_calls(bulk_shared, sc_shared, FP, 1:100, 101:200)
if(show_images) {
  show(p)
} else {
  ggsave("08_DE_FP.png", p, units = "in", dpi = 100, height = 3, width = 6)
}
p <- visualize_calls(bulk_shared, sc_shared, FN, 1:100, 101:200)
if(show_images) {
  show(p)
} else {
  ggsave("09_DE_FN.png", p, units = "in", dpi = 100, height = 3, width = 6)
}












