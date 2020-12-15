library(codaDE)
library(edgeR)

p <- 20000
sequencing_depth <- 10000
sim <- simulate_singlecell_RNAseq(p = p, n = 500, k = 1, sequencing_depth = sequencing_depth,
                                  proportion_da = 0.6, size_factor_correlation = 0, spike_in = FALSE,
                                  possible_fold_changes = NULL)

# Stringent filtering: keep things that are non-zero in at least 10% of samples
keep <- apply(sim$observed_counts, 2, function(x) {
  (sum(x > 0) / length(x)) >= 0.1
})

dgList <- DGEList(counts = t(sim$observed_counts), group = sim$groups, genes = paste0("gene", 1:p))
dgList <- dgList[keep,]

cat("No. genes profiled:",nrow(dgList),"\n")

# Re-compute the library sizes
dgList$samples$lib.size <- colSums(dgList$counts)
dgList <- calcNormFactors(dgList)
dgList <- estimateCommonDisp(dgList)
dgList <- estimateTagwiseDisp(dgList, trend = "none")

et <- exactTest(dgList)
de <- decideTestsDGE(et, adjust.method = "BH", p.value = 0.01, lfc = log(2))
summary(de)
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
