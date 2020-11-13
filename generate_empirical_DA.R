# This data is gene read counts: https://gtexportal.org/home/datasets

data_dir <- file.path("data", "GTEx_data")
data <- readRDS(file.path(data_dir, "parsed_GTEx.rds"))
annotations <- read.table(file.path(data_dir, "GTEx_annotations.txt"), header = TRUE, sep = "\t")

# Data is ~50K genes x ~17K samples
# Column 1 is Ensembl gene name, column 2 is human-readable gene name, columns 3-END are sample labels

# Remove rows with totally absent genes
rs <- rowSums(data[,3:ncol(data)])
data <- data[rs != 0,]

cat("Non-zero gene count:",nrow(data),"\n")

# Replace inscrutable sample identifiers with tissue type (from annotations)
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
  # Select a random gene
  gene_idx <- sample(1:nrow(data))[1]

  # Select two random tissue types
  ttypes <- sample(tissue_type, replace = FALSE)[1:2]

  mean_expr1 <- max(1, mean(unlist(data[gene_idx,which(colnames(data) == ttypes[1])])))
  mean_expr2 <- max(1, mean(unlist(data[gene_idx,which(colnames(data) == ttypes[2])])))
  results <- rbind(results,
                   data.frame(mean1 = mean_expr1, mean2 = mean_expr2, tissue1 = ttypes[1], tissue2 = ttypes[2]))
}

# What does the (binned) distribution over mean expression look like?
# mean_expr_distro <- c(results$mean1, results$mean2)
# table(round(mean_expr_distro, -2))

# Grab genes in percentiles of expression
# Basically half of genes have mean zero expression, so we'll bin the bottom half together
probs <- c(0, seq(from = 0.5, to = 1, length.out = 6))
boundaries <- quantile(c(results$mean1, results$mean2), probs = probs)
sub_results <- list()
for(i in 2:length(boundaries)) {
  sub_results[[i-1]] <- results[results$mean1 >= boundaries[i-1] & results$mean1 <= boundaries[i],]
}

# Save a file full of these sampled fold (quantile) changes
empirical_DA <- list()
for(i in 1:length(sub_results)) {
  empirical_DA[[i]] <- sub_results[[i]]$mean2 / sub_results[[i]]$mean1
}
saveRDS(empirical_DA, file = file.path("data", "GTEx_empirical_DA.rds"))
