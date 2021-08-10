source("path_fix.R")

library(codaDE)
library(tidyverse)
library(gridExtra)
library(randomForest)
library(optparse)

option_list = list(
  make_option(c("--dataset"),
              type = "character",
              default = "Barlow",
              help = "data set to use: VieiraSilva, Barlow, Song, Monaco, Muraro, Hashimshony, Kimmerling",
              metavar = "character"),
  make_option(c("--baseline"),
              type = "character",
              default = "self",
              help = "calls to use as a reference: self, threshold",
              metavar = "character"),
  make_option(c("--range"),
              type = "logical",
              default = "FALSE",
              help = "use bootstrapped predictions to generate an interval",
              metavar = "logical")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

dataset_name <- opt$dataset
use_baseline <- opt$baseline
use_range <- opt$range
testing <- FALSE

if(!(dataset_name %in% c("VieiraSilva", "Barlow", "Song", "Monaco", 
                         "Muraro", "Hashimshony", "Kimmerling"))) {
  stop(paste0("Invalid data set: ", dataset_name, "!\n"))
}

if(!(use_baseline %in% c("self", "threshold"))) {
  stop(paste0("Invalid DA baseline: ", use_baseline, "!\n"))
}

palette <- list(ALDEx2 = "#46A06B",
                DESeq2 = "#FF5733",
                MAST = "#EF82BB",
                NBGLM = "#7E54DE",
                scran = "#E3C012",
                simulated = "#DDDDDD")

# ------------------------------------------------------------------------------
#   Parse and wrangle validation data
# ------------------------------------------------------------------------------

if(dataset_name == "VieiraSilva") {
  abs_data <- parse_VieiraSilva(absolute = TRUE)
  rel_data <- parse_VieiraSilva(absolute = FALSE)
}
if(dataset_name == "Barlow") {
  abs_data <- parse_Barlow(absolute = TRUE)
  rel_data <- parse_Barlow(absolute = FALSE)
}
if(dataset_name == "Song") {
  abs_data <- parse_Song(absolute = TRUE)
  rel_data <- parse_Song(absolute = FALSE)
}
if(dataset_name == "Monaco") {
  abs_data <- parse_Monaco(absolute = TRUE)
  rel_data <- parse_Monaco(absolute = FALSE)
}
if(dataset_name == "Muraro") {
  abs_data <- parse_Muraro(absolute = TRUE)
  rel_data <- parse_Muraro(absolute = FALSE)
}
if(dataset_name == "Hashimshony") {
  abs_data <- parse_Hashimshony(absolute = TRUE)
  rel_data <- parse_Hashimshony(absolute = FALSE)
}
if(dataset_name == "Kimmerling") {
  abs_data <- parse_Kimmerling(absolute = TRUE)
  rel_data <- parse_Kimmerling(absolute = FALSE)
}

if(testing & nrow(abs_data$counts) > 500) {
  k <- 500
  sample_idx <- sample(1:nrow(abs_data$counts), size = k, replace = FALSE)
  abs_data$counts <- abs_data$counts[sample_idx,]
  rel_data$counts <- rel_data$counts[sample_idx,]
}

# Reorient as (samples x features)
ref_data <- t(abs_data$counts)
data <- t(rel_data$counts)
groups <- abs_data$groups

if(dataset_name == "VieiraSilva") {
  subgroups <- c("mHC", "CD")
  ref_data <- ref_data[groups %in% subgroups,]
  data <- data[groups %in% subgroups,]
  groups <- factor(groups[groups %in% subgroups])
}
if(dataset_name == "Monaco") {
  subgroups <- c("Plasmablasts", "Neutrophils")
  ref_data <- ref_data[groups %in% subgroups,]
  data <- data[groups %in% subgroups,]
  groups <- factor(groups[groups %in% subgroups])
}
if(dataset_name == "Hashimshony") {
  subgroups <- c("0", "1")
  ref_data <- ref_data[groups %in% subgroups,]
  data <- data[groups %in% subgroups,]
  groups <- factor(groups[groups %in% subgroups])
}

# Convert to integers, just for DESeq2
ref_data <- apply(ref_data, c(1,2), as.integer)
data <- apply(data, c(1,2), as.integer)

# Look at sparsity; filter out too-low features
# retain_features <- colSums(ref_data) > 10 & colSums(data) > 10
retain_features <- colMeans(ref_data) >= 1 & colMeans(data) >= 1
ref_data <- ref_data[,retain_features]
data <- data[,retain_features]

cat(paste0("Dataset dimensions: ", nrow(ref_data), " x ", ncol(ref_data), "\n"))
cat(paste0("Percent zeros: ", round(sum(data == 0)/(nrow(data)*ncol(data)), 3)*100, "%\n"))

# ------------------------------------------------------------------------------
#   Wrangle data for prediction-making
# ------------------------------------------------------------------------------

# Get info we need from data to make a prediction
if(dataset_name == "VieiraSilva") {
  counts_A <- data[groups == "mHC",]
  counts_B <- data[groups == "CD",]
}
if(dataset_name == "Barlow") {
  counts_A <- data[groups == "control",]
  counts_B <- data[groups == "keto",]
}
if(dataset_name == "Song") {
  counts_A <- data[groups == "lung",]
  counts_B <- data[groups == "brain",]
}
if(dataset_name == "Monaco") {
  counts_A <- data[groups == "Plasmablasts",]
  counts_B <- data[groups == "Neutrophils",]
}
if(dataset_name == "Muraro") {
  counts_A <- data[groups == "alpha",]
  counts_B <- data[groups == "beta",]
}
if(dataset_name == "Hashimshony") {
  counts_A <- data[groups == "0",]
  counts_B <- data[groups == "1",]
}
if(dataset_name == "Kimmerling") {
  counts_A <- data[groups == "low_mass",]
  counts_B <- data[groups == "high_mass",]
}

# This takes 2-3 min. to run on 15K features
features <- characterize_dataset(counts_A, counts_B)

# Add a few more
features$P <- ncol(counts_A)

# ------------------------------------------------------------------------------
#   Make predictions on simulations and real data and visualize these together
# ------------------------------------------------------------------------------

plot_labels <- list(FPR = "specificity (1 - FPR)", TPR = "sensitivity (TPR)")

# For the Kimmerling data, scran throws this error:
# "inter-cluster rescaling factors are not strictly positive"

#for(use_result_type in c("TPR", "FPR")) {
for(use_result_type in c("FPR")) {

  plot_df <- NULL

  # Load predictive model
  if(use_range) {
    model_fn <- file.path("output",
                          "predictive_fits",
                          "all",
                          paste0("all_", use_baseline, "_combined_", use_result_type, ".rds"))
    if(file.exists(model_fn)) {
      fit_obj_list <- readRDS(model_fn)
    } else {
      fit_files <- list.files(file.path("output", "predictive_fits", "all"),
                              pattern = paste0("all_self_.*?_", use_result_type, "\\.rds"),
                              full.names = TRUE)
      if(length(fit_files) == 0) {
        stop(paste0("Predictive model fits not found!\n"))
      }
      fit_obj_list <- list()
      for(i in 1:length(fit_files)) {
        cat(paste0("Loading predictive model ", i, " / ", length(fit_files), "\n"))
        fit_fn <- fit_files[i]
        fit_obj_list[[i]] <- readRDS(fit_fn)$result
      }
      saveRDS(fit_obj_list, model_fn)
    }
  } else {
    model_fn <- file.path("output",
                          "predictive_fits",
                          "all",
                          paste0("all_", use_result_type, "_", use_baseline, ".rds"))
    if(!file.exists(model_fn)) {
      stop(paste0("Predictive model fit not found: ", model_fn, "!\n"))
    }
    fit_obj <- readRDS(model_fn)
  }
  
  for(DE_method in c("ALDEx2", "DESeq2", "MAST", "scran")) {
    
    cat(paste0("Evaluating ", use_result_type, " x ", DE_method, "\n"))

    features$METHOD <- DE_method
    features$METHOD <- factor(features$METHOD, levels = c("ALDEx2", "DESeq2", "MAST", "scran"))

    if((max(table(groups)) < 5 | dataset_name == "Kimmerling") && DE_method == "scran") {
      next
    }
    
    # --------------------------------------------------------------------------
    #   Call discrepancy by chosen method
    # --------------------------------------------------------------------------
    
    # Get baseline differential abundance calls
    if(use_baseline == "oracle") {
      oracle_calls <- call_DA_NB(ref_data, groups)$pval
    } else {
      oracle_calls <- NULL
    }
    
    all_calls <- DA_wrapper(ref_data, data, groups, DE_method, oracle_calls)

    rates <- calc_DA_discrepancy(all_calls$calls, all_calls$oracle_calls)
    
    # --------------------------------------------------------------------------
    #   Finish feature wrangling
    # --------------------------------------------------------------------------

    # Convert to data.frame
    features_df <- as.data.frame(features)
    
    # Remove some features we're no longer using (too correlated with others)
    features_df <- features_df %>%
      select(-c(FW_RA_PFC1_D, FW_CLR_MED_D, FW_CLR_SD_D, FW_CLR_PNEG_D))
    
    # --------------------------------------------------------------------------
    #   Make predictions on simulated and real
    # --------------------------------------------------------------------------
    
    if(use_range) {
      for(fit_obj in fit_obj_list) {
        pred_real <- predict(fit_obj, newdata = features_df)
        
        plot_df <- rbind(plot_df,
                         data.frame(true = ifelse(use_result_type == "TPR", rates$TPR, 1 - rates$FPR),
                                    predicted = pred_real,
                                    type = DE_method))
      }
    } else {
      pred_real <- predict(fit_obj$result, newdata = features_df)
      
      plot_df <- rbind(plot_df,
                       data.frame(true = ifelse(use_result_type == "TPR", rates$TPR, 1 - rates$FPR),
                                  predicted = pred_real,
                                  type = DE_method))
    }
  }

  if(use_range) {
    pl <- ggplot() +
      geom_boxplot(data = plot_df,
                   mapping = aes(x = factor(type), y = predicted)) +
      geom_point(data = plot_df,
                 mapping = aes(x = factor(type), y = true, fill = type),
                 shape = 21,
                 size = 3) +
      scale_fill_manual(values = palette) +
      ylim(c(0,1)) +
      labs(x = paste0("observed ", plot_labels[[use_result_type]]),
           y = paste0("predicted ", plot_labels[[use_result_type]]),
           fill = "Data type") +
      theme(legend.position = "none")
    show(pl)
    ggsave(file.path("output",
                     "images",
                     paste0("range-validations_",
                            use_result_type,
                            "_",
                            use_baseline,
                            "-",
                            dataset_name,
                            ".png")),
           plot = pl,
           dpi = 100,
           units = "in",
           height = 4,
           width = 4)
  } else {
    pl <- ggplot() +
      geom_segment(data = data.frame(x = 0, xend = 1, y = 0, yend = 1),
                   mapping = aes(x = x, xend = xend, y = y, yend = yend)) +
      geom_point(data = plot_df,
                 mapping = aes(x = true, y = predicted, fill = type),
                 shape = 21,
                 size = 3) +
      scale_fill_manual(values = palette) + 
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      labs(x = paste0("observed ", plot_labels[[use_result_type]]),
           y = paste0("predicted ", plot_labels[[use_result_type]]),
           fill = "Data type")
    # show(pl)
    ggsave(file.path("output",
                     "images",
                     paste0("validations_",
                            use_result_type,
                            "_",
                            use_baseline,
                            "-",
                            dataset_name,
                            ".png")),
           plot = pl,
           dpi = 100,
           units = "in",
           height = 4,
           width = 5.5)
  }
}
