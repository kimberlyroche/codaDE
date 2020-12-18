# This file uses the GTEx data set to pull a big empirical sample of means and dispersions of CPM genes in a variety of tissues.
# Output of this is a bunch of files `empirical_mu_size_###.rds` that need to be assembled into one big matrix by build_baseline_model_2.R

full_data <- readRDS("empirical_mu_size_1.rds")
for(i in 2:49) {
  data <- readRDS(paste0("empirical_mu_size_",i,".rds"))
  full_data <- rbind(full_data, data)
}
saveRDS(full_data, file = "empirical_mu_size.rds")
