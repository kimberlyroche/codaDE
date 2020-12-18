# This file uses the GTEx data set to pull a big empirical sample of means and dispersions of CPM genes in a variety of tissues.
# Output of this is a bunch of files `empirical_mu_size_###.rds` that need to be assembled into one big matrix by build_baseline_model_2.R

library(MASS)
library(tidyverse)

# Function to estimate mean and dispersion of a count vector assuming an NB model
est_mu_size <- function(x) {
  fit <- tryCatch({
      glm.nb(y ~ 1, data.frame(y = x))
    },
    error = function(cond) { list() },
    warning = function(cond) { list() },
    finally = {}
  )
  if(length(fit) > 0) {
    est_mu <- unname(exp(fit$coefficients[1]))
    est_size <- fit$theta
    return(c(est_mu, est_size))
  } else {
    return(c(mean(x), NA))
  }
}

args = commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
  stop("Missing argument: tissue index!")
}

tissue_idx <- as.numeric(args[1])
testing <- TRUE

# Pull GTEx data
GTEx <- readRDS("/data/mukherjeelab/roche/codaDE/data/GTEx_data/parsed_GTEx.rds")
GTEx_annot <- read.table("/data/mukherjeelab/roche/codaDE/data/GTEx_data/GTEx_annotations.txt",
  header = TRUE, sep = "\t")
protein_coding <- read.table("mart_export.txt", header = TRUE, stringsAsFactors = FALSE)$GeneStableID

# Collect by tissue type
samples_by_tissue <- list()
for(tissue in unique(GTEx_annot$SMTSD)) {
  samples_by_tissue[[tissue]] <- GTEx_annot$SAMPID[GTEx_annot$SMTSD == tissue]
}

# Remove tissues with fewer than 50 samples
remove_idx <- c()
for(i in 1:length(samples_by_tissue)) {
  if(length(samples_by_tissue[[i]]) < 50) {
    remove_idx <- c(remove_idx, i)
  }
}
samples_by_tissue <- samples_by_tissue[-remove_idx]

# Remove `Cells - Leukemia cell line (CML)`
samples_by_tissue <- samples_by_tissue[-which(names(samples_by_tissue) == "Cells - Leukemia cell line (CML)")]
str(samples_by_tissue, max.level = 1)

# For each tissue
# (1) Pull gene expression
# (2) Collapse versions to genes
# (3) Convert to CPM
# (4) Fit a NB model to each gene and save mean and dispersion estimates

cat("Pulling versions expression data for",names(samples_by_tissue)[tissue_idx],"...\n")
gene_expr <- as.data.frame(GTEx[,which(colnames(GTEx) %in% samples_by_tissue[[tissue_idx]])])
rownames(gene_expr) <- paste(unname(unlist(GTEx[,2])), unname(unlist(GTEx[,1])), sep = "-")

# Pull protein-coding genes only by matching Ensembl IDs
gene_ids <- sapply(rownames(gene_expr), function(x) {
  ensembl_id_idx <- str_locate(x, "ENSG(\\d{11})")
  if(!is.na(ensembl_id_idx[1,1])) {
    substr(x, ensembl_id_idx[nrow(ensembl_id_idx),"start"], ensembl_id_idx[nrow(ensembl_id_idx),"end"])
  } else {
    NA
  }
})
keep_idx <- which(gene_ids %in% protein_coding)

gene_expr <- gene_expr[keep_idx,]
gene_ids <- gene_ids[keep_idx]

if(testing) {
  gene_expr <- gene_expr[sort(sample(1:nrow(gene_expr))[1:2000]),]
  gene_ids <- gene_ids[1:2000]
}

# Add Ensembl gene ID to data.frame
gene_expr_modified <- cbind(gene_label = gene_ids, gene_expr)
rownames(gene_expr_modified) <- NULL

cat("Summarizing gene expression (collapsing versions of the same gene)...\n")
# There are few duplicates; we probably don't need to do this
gene_expr_summarized <- gene_expr_modified %>%
  pivot_longer(!gene_label, names_to = "sample_id", values_to = "count") %>%
  group_by(gene_label, sample_id) %>%
  summarise(total_counts = sum(count)) %>%
  pivot_wider(names_from = sample_id, values_from = total_counts)
gene_expr_summarized <- as.data.frame(gene_expr_summarized)
rownames(gene_expr_summarized) <- gene_expr_summarized$gene_label
gene_expr_summarized <- gene_expr_summarized[,2:ncol(gene_expr_summarized)]

# Convert to CPM
gene_CPM <- apply(gene_expr_summarized, 2, function(x) {
  1e6 * (x/sum(x))
})

cat("Estimating per-gene mean and dispersion...\n")
# This will omit things that are really rare, so we're biasing our sample in that way FYI
n_genes <- nrow(gene_CPM)
if(testing) {
  n_genes <- 20
}
empirical_mu_size_tissue <- matrix(NA, n_genes, 2)
for(j in 1:n_genes) {
  x <- round(unname(gene_CPM[j,]))
  empirical_mu_size_tissue[j,] <- est_mu_size(x)
}
rownames(empirical_mu_size_tissue) <- gene_ids[1:n_genes]

saveRDS(empirical_mu_size_tissue, file = paste0("empirical_mu_size_",tissue_idx,".rds"))
