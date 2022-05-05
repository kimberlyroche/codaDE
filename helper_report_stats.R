source("path_fix.R")

datasets <- c("Barlow", "Gruen", "Hagai", "Hashimshony", "Kimmerling", "Klein",
              "Monaco", "Muraro", "Owens", "Song", "VieiraSilva", "Yu")

thresholds <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

print_names <- c("Barlow et al. (2020)",
                 'Gr\\"un et al. (2014)',
                 "Hagai et al. (2018)",
                 "Hashimshony et al. (2016)",
                 "Kimmerling et al. (2018)",
                 "Klein et al. (2015)",
                 "Monaco et al. (2019)",
                 "Muraro et al. (2016)",
                 "Owens et al. (2016)",
                 "Song et al. (2021)",
                 "Vieira-Silva et al. (2019)",
                 "Yu et al. (2014)")

blurbs <- c("single cell sequencing of zebrafish embryos; early vs. late time course samples drawn",
            "single cell RNA-sequencing of quiescent and cycling mouse fibroblasts",
            "nCounter array of human primary lung cancer vs. brain metastases",
            "single cell expression profiling of rat brain and liver tissue",
            "mouse embryonic stem cells cultured in serum and a two-inhibitor solution",
            "16S metagenomics from human gut samples of control and Crohn’s disease patients",
            "bulk RNA sequencing of both unstimulated and mock-viral infected mouse fibroblasts",
            "single cell RNA-seq of pancreatic islet cells",
            "immune cell profiling in human humans via bulk RNA-seq",
            "single cell RNA-sequencing of normally developing and leukemia inhibitory factor-treated mouse ESCs",
            "cycling, stimulated CD8+ T cells",
            "16S metagenomics from ketogenic diet and control mice")

canonical_ordering <- c(3,4,10,7,11,1,2,8,5,12,9,6)

out <- NULL

for(i in 1:length(datasets)) {
  dataset_name <- datasets[i]
  threshold <- thresholds[i]
  
  cat(paste0("Parsing data set ", print_name[i], "\n"))
  
  parsed_obj <- wrangle_validation_data(dataset_name = dataset_name,
                                        threshold = threshold,
                                        use_cpm = FALSE,
                                        testing = FALSE,
                                        hkg_list = NULL)
  ref_data <- parsed_obj$ref_data
  data <- parsed_obj$data
  groups <- parsed_obj$groups
  tax <- parsed_obj$tax
  
  counts_A <- data[groups == levels(groups)[1],]
  counts_B <- data[groups == levels(groups)[2],]
  counts_A_abs <- ref_data[groups == levels(groups)[1],]
  counts_B_abs <- ref_data[groups == levels(groups)[2],]
  
  sA <- mean(rowSums(counts_A_abs))
  sB <- mean(rowSums(counts_B_abs))
  fc <- max(sA, sB) / min(sA, sB)
  
  save_fn <- file.path("output",
                       "real_data_calls",
                       "no_norm",
                       paste0("calls_oracle_ALDEx2_",dataset_name,"_threshold",thresholds[i],"_noHKG.rds"))
  oracle_calls <- p.adjust(readRDS(save_fn)$all_calls$oracle_calls$pval, method = "BH") < 0.05
  percent_DE <- sum(oracle_calls) / length(oracle_calls)
  
  out <- rbind(out,
               data.frame(name = print_names[i],
                          blurb = blurbs[i],
                          n_feat = ncol(counts_A),
                          n_samp = paste0(sum(groups == levels(groups)[1]), ", ",
                                          sum(groups == levels(groups)[2])),
                          percent_zeros = paste0(round(sum(data == 0)/(nrow(data)*ncol(data))*100), "\\%"),
                          fc = round(fc, 1),
                          pde = paste0(round(percent_DE*100), "\\%")))
}

out <- out %>%
  arrange(pde, fc)
str_out <- ""
for(i in 1:nrow(out)) {
  str_out <- paste0(str_out, paste0(out[i,], collapse = " & "), " \\\\ \\hline \n")
}
writeLines(str_out, file.path("output", "scratch.txt"))

