library(MAST)
library(data.table)

# Using interoperability with scater given here:
# https://bioconductor.org/packages/devel/bioc/vignettes/MAST/inst/doc/MAST-interoperability.html

# And general DE method here:
# https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html#401_Visualization_of_50_most_differentially_expressed_genes

source("exploratory_analyses/other_DE/run_first.R")

cell_metadata <- data.frame(condition = as.factor(labels))
levels(cell_metadata$condition) <- c("untreated", "treated")

rownames(cell_metadata) <- colnames(count_table)
rownames(count_table) <- paste0("gene", 1:n_genes)
colnames(count_table) <- paste0("cell", 1:(n_samples_condition*2))
sce <- SingleCellExperiment(assays = list(counts = count_table),
                            colData = cell_metadata)
sce <- logNormCounts(sce)
sca = SceToSingleCellAssay(sce)

# GLM version 1
zlmCond <- zlm(~ condition, sca = sca, exprs_value = 'logcounts')

# GLM version 2 -- model out number of genes detected
# This probably shouldn't be a factor in our simulations?
# cngeneson <- colSums(assay(sca) > 0)
# colData(sca)$cngeneson <- scale(cngeneson)
# zlmCond <- zlm(~ condition + cngeneson, sca, exprs_value = 'logcounts')

summaryCond <- summary(zlmCond, doLRT = 'conditiontreated') 
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='conditiontreated' & component=='H',.(primerid, `Pr(>Chisq)`)], # hurdle P values
                  summaryDt[contrast=='conditiontreated' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') # logFC coefficients
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
head(fcHurdle)

FCTHRESHOLD <- log2(1.5)
fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)
head(fcHurdleSig)

sig_genes <- fcHurdleSig$primerid
quick_dirty_rate_calcs(sig_genes, rownames(count_table)[sim$da_genes])
