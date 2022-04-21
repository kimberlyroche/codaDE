library(codaDE)
library(phyloseq)
library(ANCOMBC)

# data <- readRDS("output/datasets/122f7064-bf7a-422e-ab4f-fbc0564684a0.rds")
data <- readRDS("output/datasets/ff4c654b-aa51-458f-9444-a5d01626e039.rds")

ref_data <- data$simulation$abundances
groups <- data$simulation$groups
data <- data$simulation$observed_counts1


sampleData <- data.frame(group = factor(groups),
                         total = rowSums(data))
# dds <- suppressMessages(DESeqDataSetFromMatrix(countData = t(data), colData = sampleData, design = ~ group + total))
dds <- suppressMessages(DESeqDataSetFromMatrix(countData = t(data), colData = sampleData, design = ~ group))
dds <- suppressMessages(DESeq2::estimateSizeFactors(object = dds))
dds <- suppressMessages(DESeq2::estimateDispersions(object = dds, fitType = "local"))
dds <- suppressMessages(DESeq2::nbinomWaldTest(object = dds))
res <- DESeq2::results(object = dds, alpha = 0.05)

pval_df <- data.frame(feature = paste0("gene", 1:ncol(data)),
                      pval = res$pvalue)
pval_df$pval[is.na(pval_df$pval)] <- 1
pval_df

plot(pval_df2$pval, pval_df$pval)

which(p.adjust(pval_df2$pval, method = "BH") < 0.05)
which(p.adjust(pval_df$pval, method = "BH") < 0.05)

plot(data[,sample(which(pval_df$pval < 0.05), size = 1)])

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
res <- dbGetQuery(conn, "SELECT FC_RELATIVE FROM datasets")
dbDisconnect(conn)

hist(res$FC_RELATIVE)




oracle_calls <- call_DA_NB(ref_data, groups)

res <- DA_wrapper(ref_data, data, groups, method = "ANCOMBC",
                  oracle_calls = oracle_calls$pval)

calc_DA_discrepancy(res$calls, res$oracle_calls)

# Define control indices
control_indices <- sample(setdiff(1:ncol(data$simulation$abundances), data$simulation$da_assignment), size = 100)

ref_data <- data$simulation$abundances
groups <- data$simulation$groups
data <- data$simulation$observed_counts1

oracle_calls <- call_DA_NB(ref_data, groups)

# Call DA with
res <- DA_wrapper(ref_data, data, groups, method = "DESeq2",
                  oracle_calls = oracle_calls$pval,
                  control_indices = control_indices)

calc_DA_discrepancy(res$calls, res$oracle_calls)

# Call DA without!
res <- DA_wrapper(ref_data, data, groups, method = "DESeq2",
                  oracle_calls = oracle_calls$pval,
                  control_indices = NULL)

calc_DA_discrepancy(res$calls, res$oracle_calls)



