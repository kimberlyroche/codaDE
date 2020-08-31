rewrite <- TRUE
metadata_file <- file.path("simulated_data", "metadata.tsv")
simdata_files <- list.files(path = "simulated_data", pattern = "*.rds")

if(!rewrite & file.exists(metadata_file)) {
  metadata <- read.table(metadata_file, header = T, stringsAsFactors = FALSE)
  keep_file_idx <- !(simdata_files %in% metadata$filename)
  simdata_files <- simdata_files[keep_file_idx]
} else {
  metadata <- data.frame(filename = c(),
    p = c(),
    proportion_da = c(),
    size_factor_correlation = c(),
    k = c(),
    filter_threshold = c(),
    in_use = c())
}

for(idx in 1:length(simdata_files)) {
  if(idx %% 1 == 0) {
    cat("Evaluating file",idx,"\n")
  }
  data <- readRDS(file.path("simulated_data", simdata_files[idx]))
  for(rr in 1:length(data$results)) {
    metadata <- rbind(metadata,
                      data.frame(filename = simdata_files[idx],
                                p = data$properties$p,
                                proportion_da = data$properties$proportion_da,
                                size_factor_correlation = data$properties$size_factor_correlation,
                                k = 1,
                                filter_threshold = data$results[[rr]]$filter_abundance,
                                in_use = FALSE)
                )
  }
}

write.table(metadata, file = metadata_file, quote = FALSE, sep = '\t', row.names = FALSE)
