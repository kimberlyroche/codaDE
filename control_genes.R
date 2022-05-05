source("path_fix.R")

library(codaDE)
library(ggplot2)
library(cowplot)

threshold <- 1
model_subdir <- "regression_plain"
use_cpm <- FALSE

out <- NULL
datasets <- c("Song", "Monaco", "Muraro", "Hagai", "Hashimshony", "Gruen")
for(dataset_name in datasets) {
  # parsed_obj <- wrangle_validation_data(dataset_name = dataset_name,
  #                                       threshold = threshold,
  #                                       use_cpm = use_cpm,
  #                                       testing = FALSE,
  #                                       hkg_list = hkg)
  # ref_data <- parsed_obj$ref_data
  # data <- parsed_obj$data
  # groups <- parsed_obj$groups
  # tax <- parsed_obj$tax
  # hkg_present <- parsed_obj$hkg_present
  # hkg_rho <- parsed_obj$hkg_rho
  
  calls <- readRDS(paste0("output/real_data_calls/no_norm/calls_oracle_DESeq2_",dataset_name,"_threshold",threshold,"_noHKG.rds"))
  calls_hkg <- readRDS(paste0("output/real_data_calls/no_norm/calls_oracle_DESeq2_CG_",dataset_name,"_threshold1.rds"))
  calls_cg <- readRDS(paste0("output/real_data_calls/no_norm/calls_oracle_DESeq2_CG_",dataset_name,"_threshold1_lowvar.rds"))
  
  out <- rbind(out,
               data.frame(name = dataset_name,
                          type = "TPR",
                          baseline = round(calls$rates$TPR, 3),
                          low_variance = round(calls_cg$rates$TPR, 3),
                          hkg = round(calls_hkg$rates$TPR, 3)))
  out <- rbind(out,
               data.frame(name = dataset_name,
                          type = "FPR",
                          baseline = round(calls$rates$FPR, 3),
                          low_variance = round(calls_cg$rates$FPR, 3),
                          hkg = round(calls_hkg$rates$FPR, 3)))
}

temp <- out %>% filter(type == "TPR")
str_out <- ""
for(i in 1:nrow(temp)) {
  str_out <- paste0(str_out,
                    temp[i,]$name,
                    " & ",
                    temp[i,]$baseline,
                    " & ",
                    temp[i,]$low_variance,
                    " & ",
                    temp[i,]$hkg,
                    " \\\\ \\hline \n")
}
writeLines(str_out, file.path("output", "scratch.txt"))

temp <- out %>% filter(type == "FPR")
str_out <- ""
for(i in 1:nrow(temp)) {
  str_out <- paste0(str_out,
                    temp[i,]$name,
                    " & ",
                    1-temp[i,]$baseline,
                    " & ",
                    1-temp[i,]$low_variance,
                    " & ",
                    1-temp[i,]$hkg,
                    " \\\\ \\hline \n")
}
writeLines(str_out, file.path("output", "scratch.txt"))











