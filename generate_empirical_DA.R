# This data is gene read counts: https://gtexportal.org/home/datasets

data_dir <- "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"
data <- readRDS(file.path(data_dir, "parsed_GTEx.rds"))
annotations <- read.table(file.path(data_dir, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"), header = TRUE, sep = "\t")

# data is ~50K genes x ~17K samples
# column 1 is Ensembl gene name, column 2 is human-readable gene name, columns 3-END are sample labels

# subset genes for testing
# subset_gene_idx <- sample(1:nrow(data))[1:5000]
# data <- data[subset_gene_idx,]
# there are some totally missing genes here; we need to remove these
rs <- rowSums(data[,3:ncol(data)])
data <- data[rs != 0,]

cat("Non-zero gene count:",nrow(data),"\n")

# replace inscrutable sample identifiers with tissue type (from annotations)
new_colnames <- unname(sapply(colnames(data)[3:ncol(data)], function(x) {
    as.character(annotations[annotations$SAMPID == x,]$SMTSD)
}))
colnames(data)[3:ncol(data)] <- new_colnames
tissue_type <- unique(new_colnames)

results <- data.frame(mean1 = c(), mean2 = c(), tissue1 = c(), tissue2 = c())
for(it in 1:100000) {
  if(it %% 5000 == 0) {
    cat("Iteration",it,"\n")
  }
  # select a random gene
  gene_idx <- sample(1:nrow(data))[1]

  # select two random tissue types
  ttypes <- sample(tissue_type, replace = FALSE)[1:2]

  mean_expr1 <- max(1, mean(unlist(data[gene_idx,which(colnames(data) == ttypes[1])])))
  mean_expr2 <- max(1, mean(unlist(data[gene_idx,which(colnames(data) == ttypes[2])])))
  results <- rbind(results,
                   data.frame(mean1 = mean_expr1, mean2 = mean_expr2, tissue1 = ttypes[1], tissue2 = ttypes[2]))
}

# what does the (binned) distribution over mean expression look like?
mean_expr_distro <- c(results$mean1, results$mean2)
table(round(mean_expr_distro, -2))

# grab genes in percentiles of expression
# basically half of genes have mean zero expression, so do graduated boundaries
probs <- c(0, seq(from = 0.5, to = 1, length.out = 6))
boundaries <- quantile(c(results$mean1, results$mean2), probs = probs)
sub_results <- list()
for(i in 2:length(boundaries)) {
  sub_results[[i-1]] <- results[results$mean1 >= boundaries[i-1] & results$mean1 <= boundaries[i],]
}

empirical_DA <- list()
for(i in 1:length(sub_results)) {
  empirical_DA[[i]] <- sub_results[[i]]$mean2 / sub_results[[i]]$mean1
}
saveRDS(empirical_DA, file = "empirical_DA.rds")

# "visualize" these distributions
#for(i in 1:length(sub_results)) {
  # before-after
  # png("condition_1-2_01.png")
  # plot(sub_results[[i]]$mean1, sub_results[[i]]$mean2, xlab = "tissue 1", ylab = "tissue 2")
  # dev.off()

#  # histogram of fold change
#  fold_change <- sub_results[[i]]$mean2 / sub_results[[i]]$mean1
#  bins <- cut(fold_change, breaks = c(0, 1/100, 1/50, 1/10, 1/4, 1/2, 1, 2, 4, 10, 50, 100, Inf))
#  print(table(bins))
#  png(paste0("condition_1-2_0",i,".png"))
#  hist(fold_change, breaks = 10)
#  dev.off()
#}
