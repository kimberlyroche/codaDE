args = commandArgs(trailingOnly=TRUE)
label <- args[1]

metadata <- read.table(file.path("simulated_data", "metadata.tsv"), header = TRUE, stringsAsFactors = FALSE)

data_dir <- file.path("simulated_analyses", label)
results_files <- list.files(path = data_dir, pattern = "*.rds")
outfile <- "results.tsv"
results <- data.frame(filename = c(), p = c(), prop_da = c(), sf_corr = c(), tp = c(), fp = c(), tn = c(), fn = c())

for(i in 1:length(results_files)) {
  if(i %% 100 == 0) {
    cat("Evaluating",i,"/",length(results_files),"\n")
  }
  data_file <- readRDS(file.path(data_dir, results_files[i]))
  results <- rbind(results,
                   data.frame(filename = results_files[i],
                              p = metadata[metadata$filename == results_files[i],]$p,
                              prop_da = metadata[metadata$filename == results_files[i],]$proportion_da,
                              sf_corr = metadata[metadata$filename == results_files[i],]$size_factor_correlation,
                              tp = data_file$tp,
                              fp = data_file$fp,
                              tn = data_file$tn,
                              fn = data_file$fn))
}

write.table(results, file = file.path(data_dir, outfile), quote = FALSE, sep = '\t', row.names = FALSE)
