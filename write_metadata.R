# query which conditions have less that 20 instances
# library(dplyr)
# m <- read.table("simulated_data/metadata.tsv", header = T, stringsAsFactors = F)
# counts <- m %>% group_by(p, proportion_da, size_factor_correlation) %>% tally()
# unfinished <- counts[which(counts$n < 20),]
# p <- as.data.frame(unfinished)
# p$n <- 20 - p$n
# for(i in 1:nrow(p)) {
#     cat(paste0("--p=",p$p[i]," --n=250 --k=1 --prop_da=",p$proportion_da[i]," --sfcorr=",p$size_factor_correlation[i],"   (",p$n[i],")\n"))
# }

# parse existing table
metadata_file <- file.path("simulated_data", "metadata.tsv")
if(file.exists(metadata_file)) {
  metadata <- read.table(metadata_file, header = T, stringsAsFactors = FALSE)
} else {
  metadata <- data.frame(filename = c(), p = c(), proportion_da = c(), size_factor_correlation = c())
}

simdata_files <- list.files(path = "simulated_data", pattern = "*.rds")
keep_file_idx <- !(simdata_files %in% metadata$filename)
simdata_files <- simdata_files[keep_file_idx]

for(idx in 1:length(simdata_files)) {
  if(idx %% 100 == 0) {
    cat("Evaluating file",idx,"\n")
  }
  data <- readRDS(file.path("simulated_data", simdata_files[idx]))
  metadata <- rbind(metadata,
                    data.frame(filename = simdata_files[idx],
                               p = data$result$p,
                               proportion_da = data$result$proportion_da,
                               size_factor_correlation = data$result$size_factor_correlation,
                               k = 1,
                               filter_threshold = 0,
                               in_use = FALSE)
              )
}

write.table(metadata, file = metadata_file, quote = FALSE, sep = '\t', row.names = FALSE)
