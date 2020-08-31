# run on 8K simulation files, this takes about 2-3 hrs.
# it's the loading of each individual data set that's rough (TBD)

output_file <- file.path("simulated_data/metadata.tsv")
metadata <- read.table(output_file, header = T, stringsAsFactors = F)

args = commandArgs(trailingOnly=TRUE)
filter_abundance <- as.numeric(args[1])
cap <- as.numeric(args[2])

metadata <- metadata[metadata$filter_threshold == filter_abundance,]

data_dir <- "simulated_data"
outfile <- paste0("results_filter_",cap,".tsv")
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
slot <- NULL
for(idx in lower_lim:upper_lim) {
  cat(paste0("Pulling data set ",idx," / ",nrow(metadata),"\n"))
  cat(paste0("\t",metadata$filename[idx],"\n"))
  data <- readRDS(file.path(data_dir, metadata$filename[idx]))
  # record the metadata we want
  if(is.null(slot)) {
    # find the results index that corresponds to this level of filtering
    slot <- 1
    while(data$results[[slot]]$filter_abundance != filter_abundance & slot <= length(data$results)) {
      slot <- slot + 1
    }
    if(slot > length(data$results)) {
      cat("Filter abundance not found in",metadata$filename[idx],"\n")
      quit()
    }
  }
  metadata[idx,"tp"] <- data$results[[slot]]$tp
  metadata[idx,"fp"] <- data$results[[slot]]$fp
  metadata[idx,"tn"] <- data$results[[slot]]$tn
  metadata[idx,"fn"] <- data$results[[slot]]$fn
  metadata[idx,"median_expr_da_quantile"] <- data$properties$median_expr_da_quantile
  metadata[idx,"net_dir_da"] <- data$properties$net_dir_da
  metadata[idx,"sparsity"] <- data$properties$sparsity
  metadata[idx,"sim_entropy"] <- data$properties$sim_entropy
  # now fix the structure of the saved object
}

metadata <- metadata[complete.cases(metadata),]

write.table(metadata, file = file.path(data_dir, outfile), quote = FALSE, sep = '\t', row.names = FALSE)
