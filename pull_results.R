args = commandArgs(trailingOnly=TRUE)
label <- args[1]
chunk_idx <- as.numeric(args[2])

metadata <- read.table(file.path("simulated_data", "metadata.tsv"), header = TRUE, stringsAsFactors = FALSE)

results_dir <- file.path("simulated_analyses", label)
data_dir <- "simulated_data"
results_files <- list.files(path = results_dir, pattern = "*.rds")
outfile <- paste0("results_",chunk_idx,".tsv")
results <- data.frame(filename = c(), p = c(), prop_da_attemped = c(), prop_da_realized = c(), sf_corr = c(), tp = c(), fp = c(), tn = c(), fn = c())

lower_idx <- max(chunk_idx, 1)
upper_idx <- min(length(results_files), chunk_idx + 199)
cat("Evaluating",lower_idx,"thru",upper_idx,"\n")
for(i in lower_idx:upper_idx) {
  if(i %% 10 == 0) {
    cat("\t",i,"/",upper_idx,"\n")
  }
  results_file <- readRDS(file.path(results_dir, results_files[i]))
  data_file <- readRDS(file.path(data_dir, results_files[i]))
  n <- round(nrow(data_file$data$observed_counts)/2)
  props_before <- rowMeans(apply(data_file$data$observed_counts[1:n,], 1, function(x) x/sum(x)))
  props_after <- rowMeans(apply(data_file$data$observed_counts[(n+1):(n*2),], 1, function(x) x/sum(x)))
  delta <- abs(props_after - props_before)
  max_prop_change <- max(delta)
  results <- rbind(results,
                   data.frame(filename = results_files[i],
                              p = metadata[metadata$filename == results_files[i],]$p,
                              prop_da_attemped = metadata[metadata$filename == results_files[i],]$proportion_da,
                              prop_da_realized = results_file$no_features_detectable / results_file$no_features_above_threshold,
                              sf_corr = metadata[metadata$filename == results_files[i],]$size_factor_correlation,
                              tp = results_file$TP,
                              fp = results_file$FP,
                              tn = results_file$TN,
                              fn = results_file$FN,
                              median_expr_da_quantile = data_file$properties$median_expr_da_quantile,
                              net_dir_da = data_file$properties$net_dir_da,
                              sparsity = data_file$properties$sparsity,
                              sim_entropy = data_file$properties$sim_entropy,
                              max_prop_change = max_prop_change))
}

write.table(results, file = file.path(results_dir, outfile), quote = FALSE, sep = '\t', row.names = FALSE)
