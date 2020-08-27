output_file <- file.path("simulated_data/metadata.tsv")
metadata <- read.table(output_file, header = T, stringsAsFactors = F)
metadata <- metadata[metadata$filter_threshold == 0,]

args = commandArgs(trailingOnly=TRUE)
cap <- as.numeric(args[1])

data_dir <- "simulated_data"
outfile <- paste0("results_nofilter_",cap,".tsv")
metadata$tp <- NULL
metadata$fp <- NULL
metadata$tn <- NULL
metadata$fn <- NULL

lower_lim <- max(1, cap - 999)
upper_lim <- min(nrow(metadata), cap)
for(idx in lower_lim:upper_lim) {
  cat(paste0("Pulling data set ",idx," / ",nrow(metadata),"\n"))
  data <- readRDS(file.path(data_dir, metadata$filename[idx]))$result
  metadata[idx,"tp"] <- data$tp
  metadata[idx,"fp"] <- data$fp
  metadata[idx,"tn"] <- data$tn
  metadata[idx,"fn"] <- data$fn
}

metadata <- metadata[complete.cases(metadata),]

write.table(metadata, file = file.path(data_dir, outfile), quote = FALSE, sep = '\t', row.names = FALSE)
