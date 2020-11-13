library(edgeR)
library(ggplot2)
library(driver)

# To do: come back and try this
# https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf

# We see compositional effects induce a lot of error once we have upwards to 1/3 - 1/2 the transcriptome
# differentially expressed.
# Use GTEx data set to determine a pseudo-max proportion of differentially expressed genes across
# different tissue types. Do we every see DE proportions this large in vivo?

GTEx <- readRDS(file.path("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct","parsed_GTEx.rds"))
GTEx_annot <- read.table(file.path("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct",
                                 "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"),
                       header = TRUE, sep = "\t")
samples_by_tissue <- list()
for(tissue in unique(GTEx_annot$SMTSD)) {
  samples_by_tissue[[tissue]] <- GTEx_annot$SAMPID[GTEx_annot$SMTSD == tissue]
}

within_tissue <- FALSE
if(within_tissue) {
  # Pull a random pair of tissues.
  tissue <- sample(names(samples_by_tissue), size = 1)
  tissue_pair <- rep(tissue, 2)
  t_idx <- which(colnames(GTEx) %in% samples_by_tissue[[tissue]])
  t1_idx <- sample(t_idx, size = 100, replace = TRUE)
  t2_idx <- sample(t_idx, size = 100, replace = TRUE)
} else {
  # Pull a random pair of tissues.
  tissue_pair <- sample(names(samples_by_tissue), size = 2, replace = FALSE)
  t1_idx <- which(colnames(GTEx) %in% samples_by_tissue[[tissue_pair[1]]])
  t2_idx <- which(colnames(GTEx) %in% samples_by_tissue[[tissue_pair[2]]])
}

cat("Using tissues:",tissue_pair[1],"x",tissue_pair[2],"\n")
gene_expr <- cbind(GTEx[,t1_idx],
                   GTEx[,t2_idx])

c_idx1 <- 1:length(t1_idx)
c_idx2 <- (length(t1_idx)+1):ncol(gene_expr)
tissues <- as.factor(c(rep("A", length(c_idx1)), rep("B", length(c_idx2))))
colnames(gene_expr) <- c(paste("A", 1:length(c_idx1), sep = "_"), paste("B", 1:length(c_idx2), sep = "_"))
rownames(gene_expr) <- paste(unname(unlist(GTEx[,2])), unname(unlist(GTEx[,1])), sep = "-")

# Filter to genes above a minimum number of non-zero appearances
keep <- which(apply(gene_expr[,c_idx1], 1, function(x) { sum(x > 0) > 5 }), apply(gene_expr[,c_idx2], 1, function(x) { sum(x > 0) > 5 }))

gene_expr <- gene_expr[keep,]

dgList <- DGEList(counts = gene_expr, genes = rownames(gene_expr))
dgList <- calcNormFactors(dgList, method = "TMM")
design <- model.matrix(~tissues)
# See also: estimateGLMCommonDisp
dgList <- estimateGLMTrendedDisp(dgList, design = design)
fit <- glmFit(dgList, design)
lrt <- glmLRT(fit, coef = 2)
is.de <- decideTestsDGE(lrt)
sum.de <- summary(is.de)
prop.de <- sum(sum.de[rownames(sum.de) != "NotSig"]) / nrow(gene_expr)
cat("Proportion DE genes:",round(prop.de, 2),"\n")
