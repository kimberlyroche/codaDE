# This file uses the GTEx data set to pull a big empirical sample of means and dispersions of CPM genes in a variety of tissues.
# Output of this is a bunch of files `empirical_mu_size_###.rds` that need to be assembled into one big matrix by build_baseline_model_2.R

file_list <- list.files(path = ".", pattern = "empirical_mu_size_pairs_\\d+_\\d+\\.rds")
full_data <- readRDS(file_list[1])
if(length(file_list) > 1) {
  for(i in 2:length(file_list)) {
    data <- readRDS(file_list[i])
    full_data <- rbind(full_data, data)
  }
}
saveRDS(full_data, file = "empirical_mu_size_pairs.rds")
