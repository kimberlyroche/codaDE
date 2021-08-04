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
              help = "data set to use: Barlow, Morton, Athanasiadou_yeast, Athanasiadou_ciona, Song, Muraro",
              metavar = "character"),
  make_option(c("--baseline"),
              type = "character",
              default = "self",
              help = "calls to use as a reference: self, threshold",
              metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

dataset_name <- opt$dataset
use_baseline <- opt$baseline

if(!(dataset_name %in% c("Barlow", "Morton", "Athanasiadou_yeast",
                         "Athanasiadou_ciona", "Song", "Muraro"))) {
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

if(dataset_name == "Barlow") {
  abs_data <- parse_Barlow(absolute = TRUE)
  rel_data <- parse_Barlow(absolute = FALSE)
}
if(dataset_name == "Morton") {
  abs_data <- parse_Morton(absolute = TRUE)
  rel_data <- parse_Morton(absolute = FALSE)
}
if(dataset_name == "Athanasiadou_ciona") {
  abs_data <- parse_Athanasiadou(absolute = TRUE, which_data = "ciona")
  rel_data <- parse_Athanasiadou(absolute = FALSE, which_data = "ciona")
  # Downsample for testing
  # k <- 1000
  # sample_idx <- sample(1:nrow(abs_data$counts), size = k, replace = FALSE)
  # abs_data$counts <- abs_data$counts[sample_idx,]
  # rel_data$counts <- rel_data$counts[sample_idx,]
}
if(dataset_name == "Athanasiadou_yeast") {
  abs_data <- parse_Athanasiadou(absolute = TRUE, which_data = "yeast")
  rel_data <- parse_Athanasiadou(absolute = FALSE, which_data = "yeast")
  # Downsample for testing
  # k <- 1000
  # sample_idx <- sample(1:nrow(abs_data$counts), size = k, replace = FALSE)
  # abs_data$counts <- abs_data$counts[sample_idx,]
  # rel_data$counts <- rel_data$counts[sample_idx,]
}
if(dataset_name == "Song") {
  # Add this to a function later if we end up using it
  file_dir <- file.path("data", "Song_2021")
  # Had to first remove a pound sign in the Accession No field name
  headers <- read.table(file.path(file_dir, "GSE161116_series_matrix.txt"),
                        header = FALSE, skip = 25, nrow = 8, sep = "\t")
  status <- unname(unlist(headers[8,2:ncol(headers)]))
  status <- factor(status, levels = c("primary lung cancer", "brain metastasis"))
  levels(status) <- c("lung", "brain")
  
  mrna <- read.table(file.path(file_dir, "GSE161116_series_matrix.txt"),
                     header = FALSE, skip = 60, nrow = 779, sep = "\t")
  rownames(mrna) <- NULL
  colnames(mrna) <- NULL
  
  gene_names <- mrna[,1]
  mrna <- mrna[,2:ncol(mrna)]
  
  ref_idx <- unname(which(sapply(gene_names, function(x) str_detect(x, "^POS_"))))
  
  # mRNA are rows 1:770
  # Negative controls (ERCC spike-ins) are named NEG_A (etc. and have very low
  #   abundance)
  # Positive controls (ERCC spike-ins) are named POS_A (etc.)
  
  # Normalize as using the spike-ins via an estimate from DESeq2
  countData <- round(mrna[ref_idx,])
  ed <- DESeqDataSetFromMatrix(countData, DataFrame(status), ~ status)
  ed <- estimateSizeFactors(ed)
  
  # Visualize the roughly-normalized data
  adjusted_mrna <- mrna[-ref_idx,]
  for(j in 1:ncol(adjusted_mrna)) {
    adjusted_mrna[,j] <- adjusted_mrna[,j] / sizeFactors(ed)[j]
  }
  
  abs_data <- list(counts = adjusted_mrna, groups = status, tax = NULL)
  rel_data <- list(counts = mrna[-ref_idx,], groups = status, tax = NULL)
}
if(dataset_name == "Muraro") {
  # Parse data and assignments
  data_orig <- read.table(file.path("data", "Muraro_2016", "GSE85241_cellsystems_dataset_4donors_updated.csv"))

  # Subset to samples with cluster assignments
  assign_fn <- file.path("data", "Muraro_2016", "cluster_assignment.rds")
  if(!file.exists(assign_fn)) {
    stop("Cluster assignments not found!")
  }
  mapping <- readRDS(assign_fn)
  data <- data_orig[,mapping %>% filter(!is.na(cluster)) %>% pull(idx)]
  # dim(data)
  
  # Pull GCG and INS sequences
  gcg_idx <- which(sapply(rownames(data), function(x) str_detect(x, "^GCG__")))
  ins_idx <- which(sapply(rownames(data), function(x) str_detect(x, "^INS__")))
  cd24_idx <- which(sapply(rownames(data), function(x) str_detect(x, "^CD24__")))
  tm4sf4_idx <- which(sapply(rownames(data), function(x) str_detect(x, "^TM4SF4__")))
  
  # ID clusters with max GCG as alpha cells, max INS as beta cells
  c1 <- data.frame(expr = unlist(unname(data[gcg_idx,])),
                   cluster = mapping$cluster[!is.na(mapping$cluster)]) %>%
    group_by(cluster) %>%
    summarize(mean_expr = mean(expr)) %>%
    arrange(desc(mean_expr)) %>%
    top_n(1) %>%
    pull(cluster)
  c2 <- data.frame(expr = unlist(unname(data[ins_idx,])), cluster = mapping$cluster[!is.na(mapping$cluster)]) %>%
    group_by(cluster) %>%
    summarize(mean_expr = mean(expr)) %>%
    arrange(desc(mean_expr)) %>%
    top_n(1) %>%
    pull(cluster)
  
  # Pull data per group
  counts_A <- data_orig[,mapping %>% filter(!is.na(cluster)) %>% filter(cluster %in% c(c1)) %>% pull(idx)]
  counts_B <- data_orig[,mapping %>% filter(!is.na(cluster)) %>% filter(cluster %in% c(c2)) %>% pull(idx)]
  counts <- cbind(counts_A, counts_B)
  
  # Define group labels
  groups <- c(rep("alpha", ncol(counts_A)), rep("beta", ncol(counts_B)))
  
  # Pull spike-in sequences and AVERAGE*
  spikein_seqs <- which(sapply(rownames(data), function(x) str_detect(x, "^ERCC-\\d+")))
  spikein_counts <- cbind(data_orig[spikein_seqs,
                                    mapping %>% filter(!is.na(cluster)) %>% filter(cluster %in% c(c1)) %>% pull(idx)],
                          data_orig[spikein_seqs,
                                    mapping %>% filter(!is.na(cluster)) %>% filter(cluster %in% c(c2)) %>% pull(idx)])
  
  # Eliminate very low count spike-ins
  spikein_counts <- spikein_counts[apply(spikein_counts, 1, function(x) min(x) >= 4),]
  size_factor <- log(unname(apply(spikein_counts, 2, mean)))
  size_factor <- scale(size_factor, scale = FALSE)
  size_factor <- size_factor * 0.25
  size_factor <- size_factor + 1
  
  # ggplot(data.frame(x = 1:length(size_factor), y = size_factor, label = groups),
  #        aes(x = x, y = y, fill = label)) +
  #   geom_point(size = 2, shape = 21)

  adj_counts <- counts
  for(i in 1:ncol(counts)) {
    adj_counts[,i] <- adj_counts[,i] / size_factor[i]
  }

  abs_data <- list(counts = adj_counts, groups = groups, tax = NULL)
  rel_data <- list(counts = counts, groups = groups, tax = NULL)
  
  # mean(colSums(abs_data$counts))
  # mean(colSums(rel_data$counts))
  
  # Downsample the observed abundances so they're in line(-ish) with the
  # reconstructed absolute abundances
  # rel_data$counts <- rel_data$counts / 30
  
  # ggplot(data.frame(x = 1:length(size_factor), y = colSums(abs_data$counts), label = groups),
  #        aes(x = x, y = y, fill = label)) +
  #   geom_point(size = 2, shape = 21)
  
  # Subset for testing
  # sample_idx <- sample(1:nrow(abs_data$counts), 500)
  # abs_data$counts <- abs_data$counts[sample_idx,]
  # rel_data$counts <- rel_data$counts[sample_idx,]
}

# Reorient as (samples x features)
ref_data <- t(abs_data$counts)
data <- t(rel_data$counts)
groups <- abs_data$groups

if(dataset_name == "Athanasiadou_ciona") {
  ref_data <- ref_data[groups %in% c("lacz", "dnfgfr"),]
  data <- data[groups %in% c("lacz", "dnfgfr"),]
  groups <- factor(groups[groups %in% c("lacz", "dnfgfr")])
}
if(dataset_name == "Athanasiadou_yeast") {
  ref_data <- ref_data[groups %in% c("C12", "C30"),]
  data <- data[groups %in% c("C12", "C30"),]
  groups <- factor(groups[groups %in% c("C12", "C30")])
}

# Convert to integers, just for DESeq2
ref_data <- apply(ref_data, c(1,2), as.integer)
data <- apply(data, c(1,2), as.integer)

# Look at sparsity; filter out too-low features
# retain_features <- colSums(ref_data) > 10 & colSums(data) > 10
retain_features <- colSums(ref_data) > 2 & colSums(data) > 2 # Loosened for Song et al.
ref_data <- ref_data[,retain_features]
data <- data[,retain_features]

cat(paste0("Dataset dimensions: ", nrow(ref_data), " x ", ncol(ref_data), "\n"))
cat(paste0("Percent zeros: ", round(sum(data == 0)/(nrow(data)*ncol(data)), 3)*100, "%\n"))

# ------------------------------------------------------------------------------
#   Wrangle data for prediction-making
# ------------------------------------------------------------------------------

# Get info we need from data to make a prediction
if(dataset_name == "Barlow") {
  counts_A <- data[groups == "control",]
  counts_B <- data[groups == "keto",]
}
if(dataset_name == "Morton") {
  counts_A <- data[groups == "before",]
  counts_B <- data[groups == "after",]
}
if(dataset_name == "Athanasiadou_ciona") {
  counts_A <- data[groups == "lacz",]
  counts_B <- data[groups == "dnfgfr",]
}
if(dataset_name == "Athanasiadou_yeast") {
  counts_A <- data[groups == "C12",]
  counts_B <- data[groups == "C30",]
}
if(dataset_name == "Song") {
  counts_A <- data[groups == "lung",]
  counts_B <- data[groups == "brain",]
}
if(dataset_name == "Muraro") {
  counts_A <- data[groups == "alpha",]
  counts_B <- data[groups == "beta",]
}

# This takes 2-3 min. to run on 15K features
features <- characterize_dataset(counts_A, counts_B)

# Add a few more
features$P <- ncol(counts_A)

# ------------------------------------------------------------------------------
#   Make predictions on simulations and real data and visualize these together
# ------------------------------------------------------------------------------

plot_labels <- list(FPR = "specificity (1 - FPR)", TPR = "sensitivity (TPR)")

for(use_result_type in c("TPR", "FPR")) {

  plot_df <- NULL
  
  for(DE_method in c("ALDEx2", "DESeq2", "MAST", "scran")) {

    features$METHOD <- DE_method
    features$METHOD <- factor(features$METHOD, levels = c("ALDEx2", "DESeq2", "MAST", "scran"))

    if(max(table(groups)) < 5 && DE_method == "scran") {
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
    #   Iterate TPR, FPR predictions
    # --------------------------------------------------------------------------
    
    model_fn <- file.path("output",
                          "predictive_fits",
                          #DE_method,
                          "all",
                          #paste0(DE_method, "_", use_result_type, "_", use_baseline, ".rds"))
                          paste0("all_", use_result_type, "_", use_baseline, ".rds"))
    
    
    if(!file.exists(model_fn)) {
      stop(paste0("Predictive model fit not found: ", model_fn, "!\n"))
    }
    fit_obj <- readRDS(model_fn) # Need to iterate DE_method too!
    
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
    
    pred_real <- predict(fit_obj$result, newdata = features_df)
    
    plot_df <- rbind(plot_df,
                     data.frame(true = ifelse(use_result_type == "TPR", rates$TPR, 1 - rates$FPR),
                                predicted = pred_real,
                                type = DE_method))
    
  }

  print(plot_df)
  
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
  show(pl)
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
