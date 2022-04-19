source("path_fix.R")

datasets <- c("VieiraSilva", "Muraro", "Hagai", "Hashimshony", "Gruen", "Kimmerling",
              "Song", "Monaco", "Yu", "Klein", "Owens", "Barlow")
thresholds <- c(rep(1, length(datasets)-1), 0)

for(i in 1:length(datasets)) {
  dataset_name <- datasets[i]
  data_obj <- readRDS(file.path("output",
                                paste0("filtered_data_",dataset_name,"_threshold",thresholds[i],".rds")))
  ref_data <- data_obj$absolute
  data <- data_obj$relative
  groups <- data_obj$groups
  if(dataset_name == "VieiraSilva") {
    counts_A <- data[groups == "mHC",]; counts_B <- data[groups == "CD",]
    counts_A_abs <- ref_data[groups == "mHC",]; counts_B_abs <- ref_data[groups == "CD",]
  }
  if(dataset_name == "Barlow") {
    counts_A <- data[groups == "control",]; counts_B <- data[groups == "keto",]
    counts_A_abs <- ref_data[groups == "control",]; counts_B_abs <- ref_data[groups == "keto",]
  }
  if(dataset_name == "Song") {
    counts_A <- data[groups == "lung",]; counts_B <- data[groups == "brain",]
    counts_A_abs <- ref_data[groups == "lung",]; counts_B_abs <- ref_data[groups == "brain",]
  }
  if(dataset_name == "Monaco") {
    counts_A <- data[groups == "CD4_naive",]; counts_B <- data[groups == "PBMC",]
    counts_A_abs <- ref_data[groups == "CD4_naive",]; counts_B_abs <- ref_data[groups == "PBMC",]
  }
  if(dataset_name == "Hagai") {
    counts_A <- data[groups == "unstimulated",]; counts_B <- data[groups == "pIC4",]
    counts_A_abs <- ref_data[groups == "unstimulated",]; counts_B_abs <- ref_data[groups == "pIC4",]
  }
  if(dataset_name == "Owens") {
    counts_A <- data[groups == "early_series",]; counts_B <- data[groups == "late_series",]
    counts_A_abs <- ref_data[groups == "early_series",]; counts_B_abs <- ref_data[groups == "late_series",]
  }
  if(dataset_name == "Klein") {
    counts_A <- data[groups == "unstimulated",]; counts_B <- data[groups == "LIF-2hr",]
    counts_A_abs <- ref_data[groups == "unstimulated",]; counts_B_abs <- ref_data[groups == "LIF-2hr",]
  }
  if(dataset_name == "Yu") {
    counts_A <- data[groups == "Brn",]; counts_B <- data[groups == "Lvr",]
    counts_A_abs <- ref_data[groups == "Brn",]; counts_B_abs <- ref_data[groups == "Lvr",]
  }
  if(dataset_name == "Muraro") {
    counts_A <- data[groups == "alpha",]; counts_B <- data[groups == "beta",]
    counts_A_abs <- ref_data[groups == "alpha",]; counts_B_abs <- ref_data[groups == "beta",]
  }
  if(dataset_name == "Hashimshony") {
    counts_A <- data[groups == "0",]; counts_B <- data[groups == "1",]
    counts_A_abs <- ref_data[groups == "0",]; counts_B_abs <- ref_data[groups == "1",]
  }
  if(dataset_name == "Kimmerling") {
    counts_A <- data[groups == "low_mass",]; counts_B <- data[groups == "high_mass",]
    counts_A_abs <- ref_data[groups == "low_mass",]; counts_B_abs <- ref_data[groups == "high_mass",]
  }
  if(dataset_name == "Gruen") {
    counts_A <- data[groups == "A",]; counts_B <- data[groups == "B",]
    counts_A_abs <- ref_data[groups == "A",]; counts_B_abs <- ref_data[groups == "B",]
  }
  
  # Report stats on this data set
  cat(paste0("Dataset: ", dataset_name, "\n"))
  cat(paste0("\tNumber of features (after filtering): ", ncol(counts_A), "\n"))
  # cat(paste0("\tNumber of samples: ", length(groups), "\n"))
  cat(paste0("\tSamples per condition (A, B): ", sum(groups == levels(groups)[1]), ", ",
             sum(groups == levels(groups)[2]), "\n"))
  cat(paste0("\tPercent zeros: ", round(sum(data == 0)/(nrow(data)*ncol(data))*100), "%\n"))
  sA <- mean(rowSums(counts_A_abs))
  sB <- mean(rowSums(counts_B_abs))
  fc <- max(sA, sB) / min(sA, sB)
  cat(paste0("\tApprox. fold change between conditions: ", round(fc, 1), "\n"))
  
  save_fn <- file.path("output",
                       "real_data_calls",
                       "no_norm",
                       paste0("calls_oracle_ALDEx2_",dataset_name,"_threshold",thresholds[i],".rds"))
  oracle_calls <- p.adjust(readRDS(save_fn)$all_calls$oracle_calls, method = "BH") < 0.05
  percent_DE <- sum(oracle_calls) / length(oracle_calls)
  cat(paste0("\tApprox. percent differential features: ", round(percent_DE*100), "%\n"))
}



