# run on 8K simulation files, this takes about 2-3 hrs.
# it's the loading of each individual data set that's rough (TBD)

output_file <- file.path("simulated_data/metadata.tsv")
metadata <- read.table(output_file, header = T, stringsAsFactors = F)
metadata <- metadata[metadata$filter_threshold == 0,]

args = commandArgs(trailingOnly=TRUE)
cap <- as.numeric(args[1])

data_dir <- "simulated_data"
outfile <- paste0("results_nofilter_",cap,".tsv")
#outfile <- "results_nofilter.tsv"

metadata$tp <- NULL
metadata$fp <- NULL
metadata$tn <- NULL
metadata$fn <- NULL
metadata$median_expr_da_quantile <- NULL
metadata$net_dir_da <- NULL
metadata$sparsity <- NULL
metadata$sim_entropy <- NULL

lower_lim <- max(1, cap - 999)
upper_lim <- min(nrow(metadata), cap)
#lower_lim <- 1
#upper_lim <- nrow(metadata)
for(idx in lower_lim:upper_lim) {
  cat(paste0("Pulling data set ",idx," / ",nrow(metadata),"\n"))
  cat(paste0("\t",metadata$filename[idx],"\n"))
  data <- readRDS(file.path(data_dir, metadata$filename[idx]))
  if("result" %in% names(data)) {
    data$properties$p <- data$result$p
    data$properties$proportion_da <- data$result$proportion_da
    data$properties$size_factor_correlation <- data$result$size_factor_correlation
    data$properties$median_expr_da_quantile <- data$result$median_expr_da_quantile
    data$properties$net_dir_da <- data$result$net_dir_da
    data$properties$sparsity <- data$result$sparsity
    data$properties$sim_entropy <- data$result$sim_entropy
    data$results <- list(data$result)
    data$result <- NULL
    data$results[[1]]$p <- NULL
    data$results[[1]]$proportion_da <- NULL
    data$results[[1]]$size_factor_correlation <- NULL
    data$results[[1]]$median_expr_da_quantile <- NULL
    data$results[[1]]$net_dir_da <- NULL
    data$results[[1]]$sparsity <- NULL
    data$results[[1]]$sim_entropy <- NULL
    saveRDS(data, file.path(data_dir, metadata$filename[idx]))
  }

  # record the metadata we want
  metadata[idx,"tp"] <- data$results[[1]]$tp
  metadata[idx,"fp"] <- data$results[[1]]$fp
  metadata[idx,"tn"] <- data$results[[1]]$tn
  metadata[idx,"fn"] <- data$results[[1]]$fn
  metadata[idx,"median_expr_da_quantile"] <- data$properties$median_expr_da_quantile
  metadata[idx,"net_dir_da"] <- data$properties$net_dir_da
  metadata[idx,"sparsity"] <- data$properties$sparsity
  metadata[idx,"sim_entropy"] <- data$properties$sim_entropy
  # now fix the structure of the saved object
}

metadata <- metadata[complete.cases(metadata),]

write.table(metadata, file = file.path(data_dir, outfile), quote = FALSE, sep = '\t', row.names = FALSE)

