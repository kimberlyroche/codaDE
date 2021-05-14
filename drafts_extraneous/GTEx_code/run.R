library(edgeR)
library(ggplot2)
library(driver)

# To do: come back and try this
# https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf

# We see compositional effects induce a lot of error once we have upwards to 1/3 - 1/2 the transcriptome
# differentially expressed.
# Use GTEx data set to determine a pseudo-max proportion of differentially expressed genes across
# different tissue types. Do we every see DE proportions this large in vivo?

# Function stolen from
slugify <- function(x, alphanum_replace="", space_replace="_", tolower=TRUE) {
  x <- gsub("[^[:alnum:] ]", alphanum_replace, x)
  x <- gsub(" ", space_replace, x)
  if(tolower) { x <- tolower(x) }
  return(x)
}

# GTEx <- readRDS("/data/mukherjeelab/roche/codaDE/data/GTEx_data/parsed_GTEx.rds")
# GTEx_annot <- read.table("/data/mukherjeelab/roche/codaDE/data/GTEx_data/GTEx_annotations.txt",
#                          header = TRUE, sep = "\t")
GTEx <- readRDS("data/GTEx_data/parsed_GTEx.rds")
GTEx_annot <- read.table("data/GTEx_data/GTEx_annotations.txt",
                         header = TRUE, sep = "\t")
samples_by_tissue <- list()
for(tissue in unique(GTEx_annot$SMTSD)) {
  samples_by_tissue[[tissue]] <- GTEx_annot$SAMPID[GTEx_annot$SMTSD == tissue]
}

# Are there differences in total abundance across tissues here?
# tissues <- names(samples_by_tissue)
# totals <- list()
# for(tissue in tissues) {
#   t_idx <- which(colnames(GTEx) %in% samples_by_tissue[[tissue]])
#   X <- GTEx[,t_idx]
#   X_totals <- unname(colSums(X))
#   totals[[tissue]] <- mean(X_totals)
# }
# 
# df <- data.frame(total = unname(unlist(totals)), tissue = names(totals))
# head(df)
# ggplot(df, aes(x = as.factor(tissue), y = total)) +
#   geom_point()

within_tissue <- FALSE
if(within_tissue) {
  # Pull a random pair of tissues.
  tissue <- sample(names(samples_by_tissue), size = 1)
  tissue_pair <- rep(tissue, 2)
  t_idx <- which(colnames(GTEx) %in% samples_by_tissue[[tissue]])
  t_all_idx <- sample(t_idx, size = min(200, length(t_idx)), replace = FALSE)
  idx <- round(length(t_all_idx)/2)
  t1_idx <- t_all_idx[1:idx]
  t2_idx <- t_all_idx[(idx+1):length(t_all_idx)]
} else {
  # Pull a random pair of tissues.
  tissue_pair <- sample(names(samples_by_tissue), size = 2, replace = FALSE)
  t1_idx <- which(colnames(GTEx) %in% samples_by_tissue[[tissue_pair[1]]])
  t2_idx <- which(colnames(GTEx) %in% samples_by_tissue[[tissue_pair[2]]])
}

pair_string <- paste0(tissue_pair[1]," x ",tissue_pair[2])
cat("Using tissues:",pair_string,"\n")
gene_expr <- cbind(GTEx[,t1_idx],
                   GTEx[,t2_idx])

c_idx1 <- 1:length(t1_idx)
c_idx2 <- (length(t1_idx)+1):ncol(gene_expr)
tissues <- as.factor(c(rep("A", length(c_idx1)), rep("B", length(c_idx2))))
colnames(gene_expr) <- c(paste("A", 1:length(c_idx1), sep = "_"), paste("B", 1:length(c_idx2), sep = "_"))
rownames(gene_expr) <- paste(unname(unlist(GTEx[,2])), unname(unlist(GTEx[,1])), sep = "-")

# Testing
# gene_expr <- gene_expr[sort(sample(1:nrow(gene_expr))[1:1000]),]

# Tutorial: https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/Supplementary-RNAseq-practical.pdf
dgList <- DGEList(counts = gene_expr, group = tissues, genes = rownames(gene_expr))
keep <- rowSums(cpm(dgList) > 5) >= 10 # CPM of at least 5 in at least 10 samples (pretty strict)
dgList <- dgList[keep,]
# Re-compute the library sizes
dgList$samples$lib.size <- colSums(dgList$counts)
dgList <- calcNormFactors(dgList)
dgList <- estimateCommonDisp(dgList)
dgList <- estimateTagwiseDisp(dgList, trend = "none")

et <- exactTest(dgList)
de <- decideTestsDGE(et, adjust.method = "BH", p.value = 0.01, lfc = log(2))
summary(de)
cat("Proportion DE genes:", round(sum(de != 0)/length(de), 2), "\n")

# Visualize the cutoff
# detags <- rownames(dgList)[as.logical(de)]
# png("test.png")
# plotSmear(et, de.tags=detags)
# abline(h=c(-1, 1), col="blue")
# dev.off()

# Visualize a results
gene <- sample(which(de == 1))[1]
png(paste0(slugify(pair_string),"_UP.png"))
plot(as.vector(unlist(dgList$counts[gene,])))
dev.off()
gene <- sample(which(de == -1))[1]
png(paste0(slugify(pair_string),"_DOWN.png"))
plot(as.vector(unlist(dgList$counts[gene,])))
dev.off()
