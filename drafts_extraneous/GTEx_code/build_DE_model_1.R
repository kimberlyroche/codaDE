# This file uses the GTEx data set to build pull an empirical sample of differences in CPM means between conditions (tissue types).

library(driver)
library(tidyverse)
library(edgeR)
library(MASS)

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

testing <- FALSE

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

# Pull a random pair of tissues and grab their bulk RNA-seq expression data
tissue_pair <- sample(names(samples_by_tissue), size = 2, replace = FALSE)
tissue_idx <- which(names(samples_by_tissue) %in% tissue_pair)
samples_by_tissue <- samples_by_tissue[tissue_pair]
str(samples_by_tissue, max.level = 1)

# Build expression matrix
t1_idx <- which(colnames(GTEx) %in% samples_by_tissue[[tissue_pair[1]]])
t2_idx <- which(colnames(GTEx) %in% samples_by_tissue[[tissue_pair[2]]])
gene_expr <- cbind(GTEx[,t1_idx],
                   GTEx[,t2_idx])
c_idx1 <- 1:length(t1_idx)
c_idx2 <- (length(t1_idx)+1):ncol(gene_expr)
tissues <- as.factor(c(rep("A", length(c_idx1)), rep("B", length(c_idx2))))
colnames(gene_expr) <- c(paste("A", 1:length(c_idx1), sep = "_"), paste("B", 1:length(c_idx2), sep = "_"))
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
  gene_expr <- gene_expr[sort(sample(1:nrow(gene_expr))[1:1000]),]
  gene_ids <- gene_ids[1:1000]
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
# gene_expr_summarized[1:5,1:5]

# Convert to CPM
gene_CPM <- apply(gene_expr_summarized, 2, function(x) {
  1e6 * (x/sum(x))
})

# Tutorial: https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/Supplementary-RNAseq-practical.pdf
dgList <- DGEList(counts = gene_CPM, group = tissues, genes = rownames(gene_CPM))
# Omit filtering
# keep <- rowSums(cpm(dgList) > 5) >= 10 # CPM of at least 5 in at least 10 samples (pretty strict)
# dgList <- dgList[keep,]
# dgList$samples$lib.size <- colSums(dgList$counts) # re-compute the library sizes
dgList <- calcNormFactors(dgList)
dgList <- estimateCommonDisp(dgList)
dgList <- estimateTagwiseDisp(dgList, trend = "none")

et <- exactTest(dgList)
de <- decideTestsDGE(et, adjust.method = "BH", p.value = 0.01, lfc = log(2))
summary(de)
cat("Proportion DE genes:", round(sum(de != 0)/length(de), 2), "\n")

# Subset and save data for just the DE-flagged genes
de_genes <- which(de != 0)
de_gene_CPM <- gene_CPM[de_genes,]
cat("Estimating per-gene mean and dispersion for DE set...\n")
# This will omit things that are really rare, so we're biasing our sample in that way FYI
n_genes <- nrow(de_gene_CPM)
if(testing) {
  n_genes <- 20
}
empirical_mu_size <- matrix(NA, n_genes, 4)
for(j in 1:n_genes) {
  x <- round(unname(gene_CPM[j,]))
  empirical_mu_size[j,1:2] <- est_mu_size(x[c_idx1])
  empirical_mu_size[j,3:4] <- est_mu_size(x[c_idx2])
}
rownames(empirical_mu_size) <- gene_ids[1:n_genes]

saveRDS(empirical_mu_size, file = paste0("empirical_mu_size_pairs_",tissue_idx[1],"_",tissue_idx[2],".rds"))
