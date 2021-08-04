source("path_fix.R")

library(tidyverse)
library(codaDE)

# ------------------------------------------------------------------------------
#   Parse and cluster Muraro et al. data as in their paper - using StemID,
#   which appears just to be a k-medoids w
# ------------------------------------------------------------------------------

# Reproduce the clustering in Muraro et al.
source("RaceID2_StemID_class.R")

data_orig <- read.table("GSE85241_cellsystems_dataset_4donors_updated.csv")

# Subset the data: 3000 features x 3072 samples (cells)
# data <- data_orig[1:3000,]
data <- data_orig

# Filter and cluster using StemID (Gruen et al.) as in Muraro et al.
sc <- SCseq(data)
sc <- filterdata(sc,
                 downsample = TRUE,
                 dsn = 1,
                 mintotal = 6000,
                 minexpr = 4,
                 minnumber = 1,
                 maxexpr = 2000,
                 rseed = 17000)
sc <- clustexp(sc,
               clustnr = 20,
               bootnr = 50,
               metric = "pearson",
               do.gap = FALSE,
               sat = TRUE,
               SE.method = "firstSEmax",
               SE.factor = .25,
               B.gap = 50,
               cln = 8,
               rseed = 17000)

# Map cluster assignments to samples
# Note: many samples seem to get kicked out in there filtering. See the
# dimensions of sc@ndata and sc@fdata versus sc@expdata
assignments <- sc@cluster$kpart
mapping <- data.frame(sample_name = colnames(data), idx = 1:ncol(data)) %>%
  left_join(data.frame(sample_name = names(assignments), cluster = unname(assignments)), by = "sample_name")
head(mapping)

saveRDS(mapping, "Muraro_cluster_assignment.rds")
quit()

# Pull spike-in sequences
# spikein_seqs <- which(sapply(rownames(data), function(x) str_detect(x, "^ERCC-\\d+")))

# Pull sequences associated with various marker genes
gcg_idx <- which(sapply(rownames(data_orig), function(x) str_detect(x, "^GCG__")))
ins_idx <- which(sapply(rownames(data_orig), function(x) str_detect(x, "^INS__")))
cd24_idx <- which(sapply(rownames(data_orig), function(x) str_detect(x, "^CD24__")))
tm4sf4_idx <- which(sapply(rownames(data_orig), function(x) str_detect(x, "^TM4SF4__")))

# Crudely differentiate cluster identity by expression of marker genes
# ggplot(data.frame(cluster = factor(assignments),
#                   expr = log2(unlist(unname(data_orig[cd24_idx,!is.na(mapping$cluster)])) + 0.5)),
#        aes(x = cluster, y = expr)) +
#   geom_boxplot()

counts_A <- data[,mapping %>% filter(!is.na(cluster)) %>% filter(cluster %in% c(5)) %>% pull(idx)]
counts_B <- data[,mapping %>% filter(!is.na(cluster)) %>% filter(cluster %in% c(7)) %>% pull(idx)]
counts <- cbind(counts_A, counts_B)
groups <- c(rep("A", ncol(counts_A)), rep("B", ncol(counts_B)))

res <- call_DA_DESeq2(t(round(counts)), groups)
x <- length(res$pval)
y <- sum(p.adjust(res$pval, method = "BH") < 0.05)
cat(paste0("Differential: ", round(y, 3), " / ", x, "\n"))













