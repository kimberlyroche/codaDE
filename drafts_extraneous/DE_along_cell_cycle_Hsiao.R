library(edgeR)

setwd("C:/Users/kim/Documents/codaDE")

# No library size normalization
call_DE_quick <- function(counts, labels) {
  cat("Proportion zeros:",sum(counts == 0)/(nrow(counts)*ncol(counts)),"\n")
  dgList <- DGEList(counts = counts, genes = rownames(counts))
  dgList <- calcNormFactors(dgList, method = "none") # alternative "TMM"
  design <- model.matrix(~labels)
  dgList <- estimateGLMTrendedDisp(dgList, design = design)
  fit <- glmFit(dgList, design)
  lrt <- glmLRT(fit, coef = 2)
  is.de <- decideTestsDGE(lrt)
  is.de
}

hsiao_data <- read.table("data/GSE121265_fucci-counts.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

lib_size <- colSums(hsiao_data[,2:ncol(hsiao_data)])
plot(lib_size, hsiao_data[1,2:ncol(hsiao_data)])

gapdh_ens <- "ENSG00000111640"
gapdh_idx <- which(hsiao_data$gene == gapdh_ens)

egfp_distro <- as.numeric(unlist(hsiao_data[1,2:ncol(hsiao_data)]))
quantiles <- quantile(egfp_distro, probs = c(0.1, 0.9))

condition1_idx <- which(egfp_distro <= quantiles[1]) + 1
condition2_idx <- which(egfp_distro >= quantiles[2]) + 1

plot(lib_size, hsiao_data[gapdh_idx,2:ncol(hsiao_data)])

counts_cond1 <- as.matrix(hsiao_data[2:nrow(hsiao_data),condition1_idx])
counts_cond2 <- as.matrix(hsiao_data[2:nrow(hsiao_data),condition2_idx])
counts <- cbind(counts_cond1, counts_cond2)
rownames(counts) <- hsiao_data$gene[2:nrow(hsiao_data)]

condition1_idx <- 1:length(condition1_idx)
condition2_idx <- (length(condition1_idx)):ncol(counts)

# Higher average expression in condition 2
mean(counts[,condition1_idx])
mean(counts[,condition2_idx])

labels <- character(ncol(counts))
labels[condition1_idx] <- "A"
labels[condition2_idx] <- "B"
labels <- as.factor(labels)

# Eliminate really lowly expressed genes
keep_genes <- unname(which(rowMeans(counts) > 1))

counts <- counts[keep_genes,]
dim(counts)

result <- call_DE_quick(counts, labels)
sum(result)/length(result)

which.de <- which(result == 1)
plot(counts[sample(which.de)[1],])

# Normalize library sizes. This should have the same effect as removing information about total abundance.
# Hard to reason about scale of discrepancy without simulating though.

counts_tpm <- apply(counts, 2, function(x) (x/sum(x))*(10^6))
result_tpm <- call_DE_quick(counts_tpm, labels)
sum(result_tpm)/length(result_tpm)

# Confusion matrix
result_v <- as.vector(result)
result_tpm_v <- as.vector(result_tpm)

# Disagreement occurs in about 20% of genes
fp <- which(result_v == 0 & result_tpm_v != 0)
fn <- which(result_v != 0 & result_tpm_v == 0)

# Is this almost entirely false negatives from normalization removing a systematic increase in abundance?
fn_samp <- sample(fn)[1]
par(mfrow = c(1,2))
plot(counts[fn_samp,])
plot(counts_tpm[fn_samp,])

