library(codaDE)
library(edgeR)

# TBD
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
  stop("Missing arguments!")
}
sequencing_depth <- as.numeric(args[1])
prop_DE <- as.numeric(args[2])

plot_gene <- function(gene, tag) {
  idx <- as.numeric(substr(gene, 5, nchar(gene)))
  keep[idx]
  # png(paste0("test_",tag,".png"), width = 800)
  par(mfrow = c(1,2))
  plot(sim$abundances[,idx])
  plot(sim$observed_counts[,idx])
  # dev.off()
}

sequencing_depth <- 50000
prop_DE <- 0.6
#p <- 20000
p <- 5000

sim <- simulate_singlecell_RNAseq(p = p, n = 500, k = 1, sequencing_depth = sequencing_depth,
                                  proportion_da = prop_DE, library_size_correlation = 0, spike_in = FALSE,
                                  possible_fold_changes = 5)

plot(rowSums(sim$abundances), rowSums(sim$observed_counts))

mean(rowSums(sim$abundances[501:1000,]))/mean(rowSums(sim$abundances[1:500,]))

# Stringent filtering: keep things that are non-zero in at least 10% of samples
keep <- apply(sim$observed_counts, 2, function(x) {
  (sum(x > 0) / length(x)) >= 0.1
})

dgList <- DGEList(counts = t(sim$observed_counts), group = sim$groups, genes = paste0("gene", 1:p))
dim(dgList)
dgList <- dgList[keep,]
dim(dgList)

cat("No. genes profiled:",nrow(dgList),"\n")

# Re-compute the library sizes
dgList$samples$lib.size <- colSums(dgList$counts)
dgList <- calcNormFactors(dgList)
dgList <- estimateCommonDisp(dgList)
dgList <- estimateTagwiseDisp(dgList, trend = "none")

et <- exactTest(dgList)
de <- decideTestsDGE(et, adjust.method = "BH", p.value = 0.01, lfc = log(2))
# summary(de)
cat("Proportion DE genes:", round(sum(de != 0)/length(de), 2), "\n")

# Of the genes that were simulated as DE, which did we catch?
# A lot of the simulated DE was in lowly abundant genes and this effectively "vanishes" with filtering, FYI.

kept_genes <- paste0("gene", which(keep))
da_sim <- paste0("gene", sim$da_genes)
da_notsim <- setdiff(paste0("gene", 1:p), paste0("gene", sim$da_genes))
da_detected <- paste0("gene", rownames(de)[de != 0])
da_notdetected <- paste0("gene", rownames(de)[de == 0])

# True positives: kept, in sim$da_genes, and non-zero in de
TP <- intersect(kept_genes, intersect(da_sim, da_detected))
# idx <- as.numeric(sapply(FN, function(x) {
#   substr(x, 5, nchar(x))
# }))

# True negatives: kept, NOT in sim$da_genes and zero in de
TN <- intersect(kept_genes, intersect(da_notsim, da_notdetected))

# False positives: kept, NOT in sim$da_genes and non-zero in de
FP <- intersect(kept_genes, intersect(da_notsim, da_detected))

# False negatives: kept, in sim$da_genes and zero in de
FN <- intersect(kept_genes, intersect(da_sim, da_notdetected))

n_TP <- length(TP)
n_TN <- length(TN)
n_FP <- length(FP)
n_FN <- length(FN)

cat("TP rate:",(n_TP / (n_TP + n_FN)),"\n")
cat("FP rate:",(n_FP / (n_FP + n_TN)),"\n")

N <- nrow(dgList)
n_TP/N
n_TN/N
n_FP/N
n_FN/N

plot_gene(sample(TP)[1], "TP")
# plot_gene(sample(TN)[1], "TN")
# plot_gene(sample(FP)[1], "FP")
# plot_gene(sample(FN)[1], "FN")




