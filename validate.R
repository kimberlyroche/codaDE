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
  make_option(c("--model"),
              type = "character",
              default = "RF",
              help = "predictive model type to use: RF, LM, or EN",
              metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

dataset_name <- opt$dataset
use_baseline <- opt$baseline
model_type <- opt$model
testing <- FALSE

if(!(dataset_name %in% c("VieiraSilva", "Barlow", "Song",
                         "Monaco", "Hagai", "Owens", "Klein", "Yu"))) {
  stop(paste0("Invalid data set: ", dataset_name, "!\n"))
}

if(!(use_baseline %in% c("self", "threshold"))) {
  stop(paste0("Invalid DA baseline: ", use_baseline, "!\n"))
}

if(!(model_type %in% c("RF", "LM", "EN"))) {
  stop(paste0("Invalid model type: ", model_type, "!\n"))
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
if(dataset_name == "Hagai") {
  abs_data <- parse_Hagai(absolute = TRUE)
  rel_data <- parse_Hagai(absolute = FALSE)
}
if(dataset_name == "Owens") {
  # This one is slow to load
  abs_data <- parse_Owens(absolute = TRUE)
  rel_data <- parse_Owens(absolute = FALSE)
}
if(dataset_name == "Klein") {
  abs_data <- parse_Klein(absolute = TRUE)
  rel_data <- parse_Klein(absolute = FALSE)
}
if(dataset_name == "Yu") {
  abs_data <- parse_Yu(absolute = TRUE)
  rel_data <- parse_Yu(absolute = FALSE)
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

# Subsample if tons of cells/samples
set.seed(100)
pairs <- table(groups)
if(pairs[1] > 100) {
  A_sample <- sample(which(groups == names(pairs)[1]), size = 100)
} else {
  A_sample <- which(groups == names(pairs)[1])
}
if(pairs[2] > 100) {
  B_sample <- sample(which(groups == names(pairs)[2]), size = 100)
} else {
  B_sample <- which(groups == names(pairs)[2])
}
ref_data <- ref_data[c(A_sample, B_sample),]
data <- data[c(A_sample, B_sample),]
groups <- groups[c(A_sample, B_sample)]

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
if(dataset_name == "Hagai") {
  counts_A <- data[groups == "unstimulated",]
  counts_B <- data[groups == "pIC4",]
}
if(dataset_name == "Owens") {
  counts_A <- data[groups == "early_series",]
  counts_B <- data[groups == "late_series",]
}
if(dataset_name == "Klein") {
  counts_A <- data[groups == "unstimulated",]
  counts_B <- data[groups == "LIF-2hr",]
}
if(dataset_name == "Yu") {
  counts_A <- data[groups == "Brn",]
  counts_B <- data[groups == "Lvr",]
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

# Pull features from existing simulations
# We'll need to scale against these later
features_sims <- pull_features(use_baseline = use_baseline,
                               exclude_partials = FALSE)
features_sims <- features_sims[sample(1:nrow(features_sims), size = 5000),]

for(use_result_type in c("TPR", "FPR")) {

  plot_df_point <- NULL
  plot_df_range <- NULL

  # Load predictive model
  model_fn <- file.path("output",
                        "predictive_fits",
                        "all",
                        paste0("all_",
                               use_baseline,
                               "_",
                               model_type,
                               "_",
                               use_result_type,
                               ".rds"))
  if(!file.exists(model_fn)) {
    stop(paste0("Predictive model fit not found: ", model_fn, "!\n"))
  }
  fit_obj <- readRDS(model_fn)

  for(DE_method in c("ALDEx2", "DESeq2", "MAST", "scran")) {
    
    cat(paste0("Evaluating ", use_result_type, " x ", DE_method, "\n"))

    features_new <- features # we'll alter this for each result type
    features_new$METHOD <- DE_method
    features_new$METHOD <- factor(features_new$METHOD,
                                  levels = c("ALDEx2", "DESeq2", "MAST", "scran"))

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
    features_new <- as.data.frame(features_new)
    
    # Remove some features we're no longer using (too correlated with others)
    features_new <- features_new %>%
      select(-c(FW_RA_PFC1_D, FW_CLR_MED_D, FW_CLR_SD_D, FW_CLR_PNEG_D))
    
    # Scale features
    # 1) Pull features from simulations
    # 2) Add this new data point
    # 3) Rescale all
    # 4) Extract this new data point
    #
    # This is tedious but I can't think of a better, simple way...

    # Remove result type we're not interested in and separate data into a
    # features and response data.frame
    features_sims2 <- features_sims %>%
      select(-one_of(ifelse(use_result_type == "TPR", "FPR", "TPR")))
    response <- features_sims2 %>%
      select(one_of(use_result_type))
    features_sims2 <- features_sims2 %>%
      select(-one_of(use_result_type))
    
    # Map in new feature
    reorder_idx <- data.frame(names = colnames(features_sims2)) %>%
      left_join(data.frame(names = colnames(features_new), idx_ds = 1:ncol(features_new)), by = "names")
    features_new <- features_new[,reorder_idx$idx_ds]
    
    features_sims2 <- rbind(features_sims2, features_new)
    
    # Scale non-factors
    factors <- which(colnames(features_sims2) %in% c("METHOD"))
    non_factors <- setdiff(1:ncol(features_sims2), factors)
    for(f in non_factors) {
      features_sims2[,f] <- scale(features_sims2[,f])
    }
    
    features_new <- features_sims2[nrow(features_sims2),]

    # --------------------------------------------------------------------------
    #   Make point predictions on simulated and real
    # --------------------------------------------------------------------------
    
    pred_real <- predict(fit_obj$result, newdata = features_new)
    
    plot_df_point <- rbind(plot_df_point,
                           data.frame(true = ifelse(use_result_type == "TPR",
                                                    rates$TPR,
                                                    1 - rates$FPR),
                                      predicted = pred_real,
                                      type = DE_method))
    
    # --------------------------------------------------------------------------
    #   Make range predictions on simulated and real
    # --------------------------------------------------------------------------

    if(model_type == "RF") {
      pred_real <- predict(fit_obj$result, newdata = features_new, predict.all = TRUE)
      
      plot_df_range <- rbind(plot_df_range,
                             data.frame(true = ifelse(use_result_type == "TPR",
                                                      rates$TPR,
                                                      1 - rates$FPR),
                                        predicted = pred_real$aggregate,
                                        type = DE_method,
                                        pred_type = "aggregate"))
      plot_df_range <- rbind(plot_df_range,
                             data.frame(true = ifelse(use_result_type == "TPR",
                                                      rates$TPR,
                                                      1 - rates$FPR),
                                        predicted = pred_real$individual[1,],
                                        type = DE_method,
                                        pred_type = "individual"))
    }
  }

  # --------------------------------------------------------------------------
  #   Visualize point predictions
  # --------------------------------------------------------------------------
  
  pl <- ggplot() +
    geom_segment(data = data.frame(x = 0, xend = 1, y = 0, yend = 1),
                 mapping = aes(x = x, xend = xend, y = y, yend = yend)) +
    geom_point(data = plot_df_point,
               mapping = aes(x = true, y = predicted, fill = type),
               shape = 21,
               size = 3) +
    scale_fill_manual(values = palette) + 
    xlim(c(0,1)) +
    ylim(c(0,1)) +
    labs(x = paste0("observed ", plot_labels[[use_result_type]]),
         y = paste0("predicted ", plot_labels[[use_result_type]]),
         fill = "Data type")
  show(pl)
  ggsave(file.path("output",
                   "images",
                   paste0("validations_",
                          model_type,
                          "_",
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
  
  # --------------------------------------------------------------------------
  #   Visualize interval predictions
  # --------------------------------------------------------------------------
  
  if(model_type == "RF") {
    pl <- ggplot() +
      geom_boxplot(data = plot_df_range[plot_df_range$pred_type == "individual",],
                   mapping = aes(x = factor(type), y = predicted),
                   width = 0.25,
                   outlier.shape = NA) +
      geom_point(data = plot_df_range[plot_df_range$pred_type == "aggregate",],
                 mapping = aes(x = factor(type), y = true, fill = type),
                 shape = 21,
                 size = 5) +
      scale_fill_manual(values = palette) +
      theme_bw() +
      ylim(c(0,1)) +
      labs(x = paste0("observed ", plot_labels[[use_result_type]]),
           y = paste0("predicted ", plot_labels[[use_result_type]]),
           fill = "Data type") +
      theme(legend.position = "none")
    show(pl)
    ggsave(file.path("output",
                     "images",
                     paste0("range-validations_",
                            model_type,
                            "_",
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
  }
}
