#' Compute reference-derived per-sample scale factor
#' 
#' @param ref_data count matrix of spike-ins (etc.) in feature x sample 
#' orientation
#' @return per-sample size factors
#' @export
compute_sf <- function(ref_data) {
  if(!is.null(dim(ref_data))) {
    # Eliminate low-count references
    # ref_data <- ref_data[rowMeans(ref_data) > 100,]
    ref_mean <- colMeans(ref_data)
  } else {
    ref_mean <- ref_data
  }
  ref_mean <- ref_mean / mean(ref_mean)
  return(ref_mean)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#
#   Vieira-Silva et al. (2019)
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Parse the Vieira-Silva et al. (2019) QMP data set
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @export
parse_VieiraSilva <- function(absolute = TRUE, use_cpm = FALSE) {
  file_dir <- file.path("data", "Vieira-Silva_2019")
  
  data <- read.table(file.path(file_dir, "QMP.matrix.tsv"))
  data <- as.matrix(data)
  lowest_observed <- min(data[data != 0])
  data <- data / lowest_observed
  
  meta <- read.table(file.path(file_dir, "41564_2019_483_MOESM3_ESM.txt"),
                     sep = "\t",
                     skip = 1,
                     header = FALSE)
  meta <- data.frame(sample_id = meta$V1, cohort = meta$V2)
  meta <- meta %>%
    filter(sample_id %in% rownames(data)) %>%
    arrange(cohort)
  meta$idx_in_meta <- 1:nrow(meta)
  
  # Get the samples in the QMP matrix into agreement
  mapping <- meta %>%
    left_join(data.frame(sample_id = rownames(data), idx_in_data = 1:nrow(data)),
              by = "sample_id")
  data <- data[mapping$idx_in_data,]
  data <- t(data)
  counts <- data
  
  groups <- factor(meta$cohort)
  
  # Subset to mHC vs. CD
  A_idx <- which(groups == "mHC")
  B_idx <- which(groups == "CD")
  counts <- counts[,c(A_idx, B_idx)]
  groups <- groups[c(A_idx, B_idx)]
  
  # Filter out ultra-low read counts
  retain_idx <- colSums(counts) > 1000
  counts <- counts[,retain_idx]
  groups <- groups[retain_idx]
  
  if(absolute) {
    if(use_cpm) {
      # Scale up these low counts
      counts <- counts * 1e06 / mean(colSums(counts))
    }
  } else {
    if(use_cpm) {
      counts <- apply(counts, 2, function(x) x/sum(x))
      counts <- round(counts*1e06)
    } else {
      # Shuffle observed abundances
      set.seed(1001)
      new_totals <- sample(round(colSums(counts)))
      for(i in 1:ncol(counts)) {
        counts[,i] <- rmultinom(1,
                                size = new_totals[i],
                                prob = counts[,i] / sum(counts[,i]))
      }
    }
  }
  
  plot(colSums(counts))
  
  counts <- data.matrix(counts)
  colnames(counts) <- NULL
  rownames(counts) <- NULL
  
  parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
  return(parsed_obj)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#
#   Barlow et al. (2020)
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Parse the Barlow et al. (2020) 16S data
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @export
parse_Barlow <- function(absolute = TRUE, use_cpm = FALSE) {
  file_dir <- file.path("data", "Barlow_2020")
  
  if(absolute) {
    data <- read.table(file.path(file_dir, "Absolute_Abundance_Table.csv"),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE)
  } else {
    data <- read.table(file.path(file_dir, "Relative_Abundance_Table.csv"),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE)
  }
    
  # Clean up data
  md_labels <- c("Diet", "Site", "Day", "mouse", "Cage")
  metadata <- data[,md_labels]
  counts <- data[,!(colnames(data) %in% md_labels)]
  counts <- counts[,2:ncol(counts)]
  
  tax <- colnames(counts)
  colnames(counts) <- NULL
  rownames(counts) <- NULL
  counts <- data.matrix(counts)
  
  # `counts` is initially 103 samples x 142 taxa
  
  # Pull stool samples from days 4, 7, 10
  keto_idx <- which(metadata$Diet == "Keto" & metadata$Day %in% c(4,7,10) & metadata$Site == "Stool")
  ctrl_idx <- which(metadata$Diet == "Control" & metadata$Day %in% c(4,7,10) & metadata$Site == "Stool")
  
  days <- metadata[c(ctrl_idx,keto_idx),]$Day
  
  counts <- rbind(counts[ctrl_idx,], counts[keto_idx,])
  labels <- c(rep("control", length(ctrl_idx)), rep("keto", length(keto_idx)))
  
  counts <- apply(counts, c(1,2), as.numeric)
  counts <- t(counts) # now features x samples
  
  # Eliminate all-zero taxa
  absent_tax_idx <- which(rowSums(counts) == 0)
  counts <- counts[-absent_tax_idx,]
  tax <- tax[-absent_tax_idx]
  
  groups <- factor(labels, levels = c("control", "keto"))
  
  if(use_cpm) {
    if(absolute) {
      # Rescale these to bring them into the same universe of scale as the 
      # relative counts
      counts <- counts * (1e6 / sum(counts[,1]))
    } else {
      # Note relative_counts are scaled to 100; scale to CPM
      counts <- counts * (1e6 / sum(counts[,1]))
    }
  }
  
  # The scaling does weird things to this data set. The minimum *observed* 
  # (non-zero) abundance is huge giving a gap between observed and unobserved 
  # features that is enormous. I'm scaling down all the counts by this minimum 
  # observed abundance for now.
  # min_observed <- min(counts[which(counts > 0, arr.ind = TRUE)])
  # counts <- counts / min_observed
  
  parsed_obj <- list(counts = counts, groups = groups, tax = tax)
  return(parsed_obj)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#
#   Song et al. (2021)
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Parse the Song et al. (2021) nCounter data
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @export
parse_Song <- function(absolute = TRUE, use_cpm = FALSE) {
  file_dir <- file.path("data", "Song_2021")
  # Had to first remove a pound sign in the Accession No field name
  headers <- read.table(file.path(file_dir, "GSE161116_series_matrix.txt"),
                        header = FALSE, skip = 25, nrow = 8, sep = "\t")
  groups <- unname(unlist(headers[8,2:ncol(headers)]))
  groups <- factor(groups, levels = c("primary lung cancer", "brain metastasis"))
  levels(groups) <- c("lung", "brain")
  
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
  
  # Observed abundances correlate very well to spike-in abundances already;
  # forego renormalization

  counts <- data.matrix(mrna)
  if(use_cpm) {
    counts <- apply(counts, 2, function(x) x/sum(x))
    counts <- round(counts*1e06)
  }

  counts <- counts[-ref_idx,]
  if(!absolute) {
    # Shuffle observed abundances
    set.seed(1000)
    new_totals <- sample(round(colSums(counts)))
    for(i in 1:ncol(counts)) {
      counts[,i] <- rmultinom(1, size = new_totals[i],
                              prob = counts[,i] / sum(counts[,i]))
    }
  }
  colnames(counts) <- NULL
  rownames(counts) <- NULL
  
  parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
  return(parsed_obj)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#
#   Monaco et al. (2019)
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Parse the Monaco et al. (2019) immune cell data set
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @import edgeR
#' @export
parse_Monaco <- function(absolute = TRUE, use_cpm = FALSE) {
  file_dir <- file.path("data", "Monaco_2019")
  
  # Parse data and assignments
  data <- read.table(file.path(file_dir,
                               "GSE107011_Processed_data_TPM.txt"))
  
  spike_idx <- which(sapply(rownames(data), function(x) str_detect(x, "^ERCC-")))
  spike_counts <- data[spike_idx,]
  
  # Eliminate some samples where spike-in counts are super low. These look like
  # fishy, low quality samples.
  spike_sums <- colSums(spike_counts)
  retain_idx <- which(spike_sums >= quantile(spike_sums, probs = c(0.05)))
  data <- data[,retain_idx]

  if(use_cpm) {
    data <- apply(data, 2, function(x) x/sum(x))
    data <- round(data*1e06)
  }

  spike_counts <- data[spike_idx,]

  # Plot spike-in relative abundance by type
  groups <- sapply(colnames(data), function(x) {
    pieces <- strsplit(x, "_")[[1]]
    paste0(pieces[2:length(pieces)], collapse = "_")
  })
  
  A_idx <- which(groups == "PBMC")
  B_idx <- which(groups == "CD4_naive")
  data <- data[,c(A_idx,B_idx)]
  spike_counts <- spike_counts[,c(A_idx,B_idx)]
  groups <- groups[c(A_idx,B_idx)]
  
  # Re-normalize
  sf <- compute_sf(data[spike_idx,])
  counts <- data[-spike_idx,]
  if(absolute) {
    for(i in 1:ncol(counts)) {
      counts[,i] <- counts[,i] / sf[i]
    }
  } else {
    # Shuffle observed abundances
    set.seed(101)
    new_totals <- sample(round(colSums(counts)))
    for(i in 1:ncol(counts)) {
      counts[,i] <- rmultinom(1, size = new_totals[i],
                              prob = counts[,i] / sum(counts[,i]))
    }
  }
  feature_names <- rownames(counts)
  counts <- data.matrix(counts)
  colnames(counts) <- NULL
  rownames(counts) <- NULL
  
  reorder <- order(groups)
  counts <- counts[,reorder]
  groups <- groups[reorder]
  groups <- unname(groups)
  
  # View differences in average totals across cell types
  # palette <- generate_highcontrast_palette(40)
  # plot_df <- data.frame(total = colSums(counts),
  #                       type = groups) %>%
  #   group_by(type) %>%
  #   mutate(mid_count = median(total)) %>%
  #   select(-total) %>%
  #   distinct()
  # ggplot(plot_df, aes(x = 1:nrow(plot_df), y = mid_count, fill = factor(type))) +
  #   geom_bar(stat = "identity") +
  #   scale_fill_manual(values = palette) +
  #   labs(x = "type index",
  #        y = "median total abundance",
  #        fill = "Cell type")
  
  # Subset to PBMCs vs. naive CD4 cells
  # A_idx <- which(groups == "PBMC")
  # B_idx <- which(groups == "CD4_naive")
  # counts <- counts[,c(A_idx, B_idx)]
  # groups <- unname(groups[c(A_idx, B_idx)])
  
  parsed_obj <- list(counts = counts, groups = groups, tax = feature_names)
  return(parsed_obj)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#
#   Hagai et al. (2018)
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Parse the Hagai et al. (2018) cross-species, cross-cell, multiple condition
#' experiment for mouse unstimulated vs. pIC4 time point
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @export
parse_Hagai <- function(absolute = TRUE, use_cpm = FALSE) {
  counts_A <- read.table(file.path("data",
                                   "Hagai_2018",
                                   "counts_mouse_unst.txt"),
                         header = TRUE)
  counts_B <- read.table(file.path("data",
                                   "Hagai_2018",
                                   "counts_mouse_pIC4.txt"),
                         header = TRUE)
  
  # Genes match between datasets
  gene_IDs <- counts_A$gene
  
  counts_A <- counts_A[,2:ncol(counts_A)]
  counts_B <- counts_B[,2:ncol(counts_B)]
  
  # Eliminate lowly sequenced samples
  # counts_A <- counts_A[,colSums(counts_A) > 5000]
  # counts_B <- counts_B[,colSums(counts_B) > 5000]
  
  counts <- round(cbind(counts_A, counts_B))

  if(use_cpm) {
    counts <- apply(counts, 2, function(x) x/sum(x))
    counts <- round(counts*1e06)
  }

  batch <- unname(sapply(colnames(counts), function(x) {
    str_split(x, "_")[[1]][1]
  }))
  condition <- c(rep("unstimulated", ncol(counts_A)),
                 rep("pIC4", ncol(counts_B)))
  
  # Pull spike-in sequences
  ref_idx <- which(sapply(gene_IDs, function(x) str_detect(x, "^ERCC-")))
  spike_counts <- counts[ref_idx,]
  sf <- compute_sf(spike_counts)
  counts <- counts[-ref_idx,]
  if(absolute) {
    for(j in 1:ncol(counts)) {
      counts[,j] <- counts[,j]/sf[j]
    }
  } else {
    # Shuffle observed abundances
    set.seed(1000)
    new_totals <- sample(round(colSums(counts)))
    for(i in 1:ncol(counts)) {
      counts[,i] <- rmultinom(1, size = new_totals[i],
                              prob = counts[,i] / sum(counts[,i]))
    }
  }
  rownames(counts) <- NULL
  colnames(counts) <- NULL
  counts <- data.matrix(counts)
  groups <- condition
  
  # totals1 <- colSums(counts)[groups == unique(groups)[1]]
  # sd1 <- sd(totals1)
  # outliers1 <- which((totals1 > mean(totals1) + sd1*2) | (totals1 < mean(totals1) - sd1*2))
  # totals2 <- colSums(counts)[groups == unique(groups)[2]]
  # sd2 <- sd(totals2)
  # outliers2 <- which((totals2 > mean(totals2) + sd2*2) | (totals2 < mean(totals2) - sd2*2))
  # outliers <- c(outliers1, outliers2)
  
  # if(length(outliers) > 0) {
  #   # TBD
  # }
  
  parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
  return(parsed_obj)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#
#   Owens et al. (2016)
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Parse the Owens et al. (2016) zebrafish developmental series
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @export
parse_Owens <- function(absolute = TRUE, use_cpm = FALSE) {
  counts <- read.table(file.path("data",
                                 "Owens_2016",
                                 "clutchA_polya_relative_TPM_gene_isoform.txt"),
                       header = FALSE,
                       fill = TRUE)
  timeline <- as.numeric(unname(unlist(counts[2,1:(ncol(counts)-2)])))
  gene_IDs <- counts[3:nrow(counts),2]
  counts <- counts[3:nrow(counts),3:ncol(counts)]
  rownames(counts) <- NULL
  colnames(counts) <- NULL
  counts <- apply(counts, c(1,2), as.numeric) # This takes ~1 min.
  
  # Pull spike-in sequences
  ref_idx <- which(sapply(gene_IDs, function(x) str_detect(x, "^ERCC-")))
  
  if(use_cpm) {
    counts <- apply(counts, 2, function(x) x/sum(x))
    counts <- round(counts*1e06)
  }
  
  spike_counts <- counts[ref_idx,]
  counts <- counts[-ref_idx,]
  sf <- compute_sf(spike_counts)

  if(absolute) {
    for(j in 1:ncol(counts)) {
      counts[,j] <- counts[,j]/sf[j]
    }
  } else {
    # Shuffle observed abundances
    set.seed(1000)
    new_totals <- sample(round(colSums(counts)))
    for(i in 1:ncol(counts)) {
      counts[,i] <- rmultinom(1, size = new_totals[i],
                              prob = counts[,i] / sum(counts[,i]))
    }
  }

  A_idx <- which(timeline < 12)
  B_idx <- which(timeline > 30)
  counts <- counts[,c(A_idx, B_idx)]
  groups <- c(rep("early_series", length(A_idx)), rep("late_series", length(B_idx)))
  rownames(counts) <- NULL
  colnames(counts) <- NULL

  parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
  return(parsed_obj)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#
#   Klein et al. (2015)
#
#   GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65525
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Parse the Klein et al. (2015) mouse +/- LIF data set
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @export
parse_Klein <- function(absolute = TRUE, use_cpm = FALSE) {
  counts_A <- read.table(file.path("data",
                                   "Klein_2015",
                                   "GSM1599494_ES_d0_main.csv",
                                   "GSM1599494_ES_d0_main.csv"),
                         sep = ",")
  counts_B <- read.table(file.path("data",
                                   "Klein_2015",
                                   "GSM1599497_ES_d2_LIFminus.csv",
                                   "GSM1599497_ES_d2_LIFminus.csv"),
                         sep = ",")
  # Gene IDs same order
  gene_IDs <- counts_A$V1
  counts_A <- counts_A[,2:ncol(counts_A)]
  counts_B <- counts_B[,2:ncol(counts_B)]
  
  gapdh_idx <- which(gene_IDs == "Gapdh")
  counts <- cbind(counts_A, counts_B)
  groups <- c(rep("unstimulated", ncol(counts_A)), rep("LIF-2hr", ncol(counts_B)))
  
  # Remove samples with low counts
  gapdh_expr <- unlist(counts[gapdh_idx,])
  remove_idx <- which(colSums(counts) < 5000 | gapdh_expr < 1)
  counts <- counts[,-remove_idx]
  groups <- groups[-remove_idx]
  
  gapdh_expr <- unlist(counts[gapdh_idx,])
  sf <- compute_sf(gapdh_expr)

  if(use_cpm) {
    counts <- apply(counts, 2, function(x) x/sum(x))
    counts <- round(counts*1e06)
  }

  counts <- counts[-gapdh_idx,]

  if(absolute) {
    for(j in 1:ncol(counts)) {
      counts[,j] <- counts[,j]/sf[j]
    }
  } else {
    # Shuffle observed abundances
    set.seed(1000)
    new_totals <- sample(round(colSums(counts)))
    for(i in 1:ncol(counts)) {
      counts[,i] <- rmultinom(1, size = new_totals[i],
                              prob = counts[,i] / sum(counts[,i]))
    }
  }
  rownames(counts) <- NULL
  colnames(counts) <- NULL
  counts <- data.matrix(counts)

  parsed_obj <- list(counts = counts, groups = groups, tax = gene_IDs[-gapdh_idx])
  return(parsed_obj)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#
#   Yu et al. (2014)
#
#   GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53960
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Parse the Yu et al. (2014) rat tissue bulk RNA-seq
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @export
parse_Yu <- function(absolute = TRUE, use_cpm = FALSE) {
  # type <- c("Adr", "Brn", "Hrt", "Kdn", "Lng", "Lvr", "Msc", "Spl", "Thm", "Tst", "Utr")
  # Brain and liver are the most differential in terms of sample abundance
  type <- c("Brn", "Lvr")
  counts <- NULL
  groups <- c()
  for(tt in type) {
    files <- list.files(file.path("data", "Yu_2014"),
                        pattern = paste0("GSM1328.*?_SEQC_",tt,"_.*?.txt$"),
                        full.names = TRUE)
    for(i in 1:length(files)) {
      groups <- c(groups, tt)
      if(is.null(counts)) {
        counts <- read.table(files[i], header = TRUE)
      } else {
        counts2 <- read.table(files[i], header = TRUE)
        counts <- cbind(counts, unlist(counts2[,2]))
      }
    }
  }
  gene_IDs <- counts$AceVeiwGeneSymbol
  rownames(counts) <- NULL
  colnames(counts) <- NULL
  counts <- counts[,2:ncol(counts)]
  counts <- apply(counts, c(1,2), as.numeric)

  # Per Supplementary Data 1 in their paper, (log) total observed abundances 
  # ("log2FPKM") appear to very closely correlate with (log) total ERCC 
  # abundances ("log2(ERCC)"). That's how I'm interpreting their data: as
  # accurate total abundances.
  # See: https://www.nature.com/articles/ncomms4230#accession-codes
  
  if(absolute) {
    if(use_cpm) {
      counts <- counts * 1e06 / mean(colSums(counts))
    }
  } else {
    if(use_cpm) {
      if(use_cpm) {
        counts <- apply(counts, 2, function(x) x/sum(x))
        counts <- round(counts*1e06)
      }
      # Shuffle observed abundances
      set.seed(1000)
      new_totals <- sample(round(colSums(counts)))
      for(i in 1:ncol(counts)) {
        counts[,i] <- rmultinom(1, size = new_totals[i],
                                prob = counts[,i] / sum(counts[,i]))
      }
    }
  }
  counts <- data.matrix(counts)

  parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
  return(parsed_obj)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#
#   Morton et al. (2019)
#
#   No longer in use!
#   Very little detectable differential abundance in "absolute" counts here
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Parse the Morton et al. (2019) 16S data
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @export
parse_Morton <- function(absolute = TRUE) {
  file_dir <- file.path("data", "Morton_2019", "Github")
  
  otu_table <- read.delim(file.path(file_dir, "oral_trimmed_deblur.txt"),
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  metadata <- read.delim(file.path(file_dir, "oral_trimmed_metadata.txt"),
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  if(absolute) {
    # "flowcount" is calculated this way in the paper
    # (See their Python notebook for Fig. 2)
    flowcount <- (metadata$flow.cell.5min.1 + metadata$flow.cell.5min.2) / 2
  }
  event_labels <- metadata$brushing_event
  subject_labels <- metadata$HostSubject
  timepoint_labels <- metadata$`Timepoint.`
  timepoint_labels[timepoint_labels %in% c(1,2)] <- "_AM"
  timepoint_labels[timepoint_labels %in% c(7,8,9)] <- "_PM"
  sample_labels <- metadata$SampleID
  
  # Compute relative abundances
  tax_legend <- otu_table$OTU_ID
  
  if(absolute) {
    props <- otu_table[,2:ncol(otu_table)]
    props <- apply(props, 2, function(x) x/sum(x))
    
    # Scale these up to absolute abundances
    # `props` is taxa x samples
    counts <- props
    for(i in 1:ncol(counts)) {
      counts[,i] <- counts[,i]*flowcount[i]
    }
  } else {
    counts <- otu_table[,2:ncol(otu_table)]
  }
  
  # Re-order samples by condition (event before/after) and subject (where most
  # subjects have a morning and evening brushing event)
  before_events <- which(event_labels == "before")
  before_subjects <- paste0(subject_labels[before_events],
                            timepoint_labels[before_events])
  after_events <- which(event_labels == "after")
  after_subjects <- paste0(subject_labels[after_events],
                           timepoint_labels[after_events])
  counts <- counts[,c(before_events, after_events)]
  subjects <- factor(c(before_subjects, after_subjects))
  
  groups <- factor(c(rep("before", length(before_events)),
                     rep("after", length(after_events))),
                   levels = c("before", "after"))
  
  if(absolute) {
    # Scale by minimum observed abundance as I did with Barlow et al.
    min_observed <- min(counts[which(counts > 0, arr.ind = TRUE)])
    counts <- counts / min_observed
  }
  
  counts <- data.matrix(counts)
  rownames(counts) <- NULL
  colnames(counts) <- NULL
  
  parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
  return(parsed_obj)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#
#   Athanasiadou et al. (2021)
#
#   No longer in use!
#   Very little fold change between conditions here
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Utility function that applies the ERCC normalization code from Athanasiadou
#' et al. (2021)
#'
#' @param file_dir base file directory
#' @param SI_file experiment SI file
#' @param RNA_file experiment RNA file
#' @param ERCC_annot experiment ERCC file
#' @param groups per-sample group labels
#' @return matrix of rescaled count data
#' @export
normalize_Athanasiadou <- function(file_dir, SI_file, RNA_file,
                                   ERCC_annot, groups) {
  # Spike-in counts and RNA counts
  YmatSI <- as.matrix(read.table(file = file.path(file_dir, SI_file)))
  YmatRNA <- as.matrix(read.table(file = file.path(file_dir, RNA_file)))
  
  # Spike-in amounts added (in attomoles)
  ERCC <- read.table(file = file.path(file_dir, ERCC_annot))
  ERCC <- ERCC[rownames(YmatSI),]
  
  # Spike-in library sizes
  LibSizeSI <- colSums(YmatSI)
  
  # RNA library sizes
  LibSizeRNA <- colSums(YmatRNA)
  
  Ratio <- LibSizeSI/LibSizeRNA
  # data.frame(LibSizeSI, LibSizeRNA, Ratio) # print
  
  # Perform the normalization as per RMarkdown file
  # Total amol of each spike-in
  Nvec <- ERCC[,"attomoles"]
  # In sample to be sequenced
  names(Nvec) <- rownames(ERCC)
  AttomoleToMoleculesPerCell <- (1.e-18)*(6.22e23)*(1.e-7)
  MoleculesPerCell <- Nvec*AttomoleToMoleculesPerCell
  names(MoleculesPerCell) <- names(Nvec)
  
  # From RMarkdown comments: "Compute spike-in proportions. Choose as the 
  # reference spike-in, the one with the largest proportion; record its 
  # proportion and attomoles; and compute the library-specific $\nu_j$ 
  # calibration factors"
  f.vec <- rowSums(YmatSI)/sum(YmatSI) # fraction of total spike-in
  # Counts accounted for by each spike-in
  INDmax <- which(f.vec == max(f.vec)) # index of spike-in
  f.ref <- f.vec[INDmax] # fraction (empirical proportion) for reference
  n.ref <- Nvec[names(f.vec[INDmax])] # corresponding amol
  nu <- (f.ref/n.ref)*LibSizeSI # for counts to nominal amol conversion
  Zmat <- YmatRNA%*%diag(1/nu) # normalization giving amol
  colnames(Zmat) <- colnames(YmatRNA)
  return(Zmat)
}

#' Parse the Athanasiadou et al. (2021) bulk RNA-seq data
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @param which_data options are "ciona" or "yeast"
#' @return named list of counts and group labels
#' @export
parse_Athanasiadou <- function(absolute = TRUE, which_data = "ciona") {
  which_data <- tolower(which_data)
  if(!(which_data %in% c("ciona", "yeast"))) {
    stop("Invalid data type!")
  }
  if(which_data == "yeast") {
    file_dir <- file.path("data", "Athanasiadou_2021", "S1CodeandData")
    
    groups <- as.factor(rep(c("C12", "C20", "C30"), rep(3,3)))
    
    if(absolute) {
      counts <- normalize_Athanasiadou(file_dir = file_dir,
                                       SI_file = "YmatSIyeast.txt",
                                       RNA_file = "YmatRNAyeast.txt",
                                       ERCC_annot = "ERCCyeast.txt",
                                       groups = groups)
    } else {
      counts <- as.matrix(read.table(file = file.path(file_dir,
                                                      RNA_file = "YmatRNAyeast.txt")))
    }
    
    parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
    # saveRDS(parsed_obj, file = parsed_data_fn)
    return(parsed_obj)
  } else {
    # Ciona robusta
    file_dir <- file.path("data", "Athanasiadou_2021", "S1CodeandData")
    
    groups <- as.factor(rep(c("lacz", "dnfgfr", "camras"), rep(3,3)))
    
    if(absolute) {
      counts <- normalize_Athanasiadou(file_dir = file_dir,
                                       SI_file = "YmatSIciona.txt",
                                       RNA_file = "YmatRNAciona.txt",
                                       ERCC_annot = "ERCCciona.txt",
                                       groups = groups)
    } else {
      counts <- as.matrix(read.table(file = file.path(file_dir,
                                                      RNA_file = "YmatRNAciona.txt")))
    }
    
    baseline_counts <- counts[,which(groups == "lacz")]
    absent_tax_idx <- which(rowSums(baseline_counts) == 0)
    baseline_counts <- baseline_counts[,-absent_tax_idx]
    mean_ciona <- mean(baseline_counts)
    
    # Ciona counts are super low relative to (abundant) spike-in
    # Scale these up to the mean counts in the S. cerevisieae experiment
    mean_yeast <- 7.921205 # I've pre-calculated this, as above
    scale <- mean_yeast / mean_ciona
    counts <- counts * scale
    
    rownames(counts) <- NULL
    colnames(counts) <- NULL
    counts <- data.matrix(counts)
    
    # Subset to the SC samples
    A_idx <- which(groups == "dnfgfr")
    B_idx <- which(groups == "lacz")
    counts <- counts[,c(A_idx, B_idx)]
    groups <- c(rep("dnfgfr", length(A_idx)), rep("lacz", length(B_idx)))
    
    parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
    return(parsed_obj)
  }
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#
#   Muraro et al. (2016)
#
#   No longer in use!
#   Very little fold change between conditions here
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Parse the Muraro et al. (2016) single-cell pancreas data set
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @import dplyr
#' @export
parse_Muraro <- function(absolute = TRUE, use_cpm = FALSE) {
  file_dir <- file.path("data", "Muraro_2016")
  
  # Parse data and assignments
  data_orig <- read.table(file.path(file_dir,
                                    "GSE85241_cellsystems_dataset_4donors_updated.csv"))
  
  # Subset to samples with cluster assignments
  # The authors of the paper don't give cell type labels in their supplemental
  # data. They uses a k-medoids clustering algorithm (StemID), then assign 
  # likely identity to clusters based on marker gene expression.
  # We'll do a quick and dirty version of the same.
  assign_fn <- file.path(file_dir, "cluster_assignment.rds")
  if(!file.exists(assign_fn)) {
    stop("Cluster assignments not found!")
  }
  mapping <- readRDS(assign_fn)
  data <- data_orig[,mapping %>% filter(!is.na(cluster)) %>% pull(idx)]
  
  # Pull GCG and INS sequences
  gcg_idx <- which(sapply(rownames(data), function(x) str_detect(x, "^GCG__")))
  ins_idx <- which(sapply(rownames(data), function(x) str_detect(x, "^INS__")))
  # cd24_idx <- which(sapply(rownames(data), function(x) str_detect(x, "^CD24__")))
  # tm4sf4_idx <- which(sapply(rownames(data), function(x) str_detect(x, "^TM4SF4__")))
  
  # ID clusters with max GCG as alpha cells, max INS as beta cells
  c1 <- suppressMessages(data.frame(expr = unlist(unname(data[gcg_idx,])),
                                    cluster = mapping$cluster[!is.na(mapping$cluster)]) %>%
                           group_by(cluster) %>%
                           summarize(mean_expr = mean(expr)) %>%
                           arrange(desc(mean_expr)) %>%
                           top_n(1) %>%
                           pull(cluster))
  c2 <- suppressMessages(data.frame(expr = unlist(unname(data[ins_idx,])),
                                    cluster = mapping$cluster[!is.na(mapping$cluster)]) %>%
                           group_by(cluster) %>%
                           summarize(mean_expr = mean(expr)) %>%
                           arrange(desc(mean_expr)) %>%
                           top_n(1) %>%
                           pull(cluster))
  
  # Pull data per group
  counts_A <- data_orig[,mapping %>% filter(!is.na(cluster)) %>% filter(cluster %in% c(c1)) %>% pull(idx)]
  counts_B <- data_orig[,mapping %>% filter(!is.na(cluster)) %>% filter(cluster %in% c(c2)) %>% pull(idx)]
  counts <- cbind(counts_A, counts_B)
  
  # Define group labels
  groups <- c(rep("alpha", ncol(counts_A)), rep("beta", ncol(counts_B)))
  
  # Pull spike-in sequences
  spikein_seqs <- which(sapply(rownames(data), function(x) str_detect(x, "^ERCC-\\d+")))
  spikein_counts <- cbind(data_orig[spikein_seqs,
                                    mapping %>% filter(!is.na(cluster)) %>% filter(cluster %in% c(c1)) %>% pull(idx)],
                          data_orig[spikein_seqs,
                                    mapping %>% filter(!is.na(cluster)) %>% filter(cluster %in% c(c2)) %>% pull(idx)])
  
  # Use edgeR (for now) to compute the size factor
  # sf <- calcNormFactors(spikein_counts) # columns assumed to be samples
  sf <- compute_sf(spikein_counts)
  
  if(use_cpm) {
    counts <- apply(counts, 2, function(x) x/sum(x))
    counts <- round(counts*1e06)
  }
  
  counts <- counts[-spikein_seqs,]
  
  # Not strong
  # plot_df <- data.frame(total = colSums(counts),
  #                       group = groups,
  #                       reference = sf)
  # plot_df <- plot_df %>%
  #   filter(total > quantile(plot_df$total, probs = c(0.05)) & total < quantile(plot_df$total, probs = c(0.95))) %>%
  #   filter(reference > quantile(plot_df$reference, probs = c(0.05)) & reference < quantile(plot_df$reference, probs = c(0.95)))
  # ggplot(plot_df, aes(x = reference, y = total)) +
  #   geom_smooth(aes(color = group), formula = y ~ x, method = "lm", alpha = 0.5) +
  #   geom_point(aes(fill = group), size = 2, shape = 21) +
  #   scale_fill_manual(values = c("#58D68D", "#9B59B6")) +
  #   scale_color_manual(values = c("#58D68D", "#9B59B6"), guide = FALSE) +
  #   theme_bw() +
  #   labs(x = "Reference (Gadph)", y = "Total abundance", fill = "Group")
  # ggsave(file.path("output", "images", "Muraro_totals.png"),
  #        units = "in",
  #        dpi = 100,
  #        height = 5,
  #        width = 7)
  # for(grp in unique(groups)) {
  #   r2 <- cor(plot_df[plot_df$group == grp,]$reference,
  #             plot_df[plot_df$group == grp,]$total)**2
  #   cat(paste0("R^2 (", grp, "): ", round(r2, 3), "\n"))
  # }
  
  if(absolute) {
    for(j in 1:ncol(counts)) {
      counts[,j] <- counts[,j] / sf[j]
    }
  } else {
    set.seed(1000)
    new_totals <- sample(round(colSums(counts)))
    for(i in 1:ncol(counts)) {
      counts[,i] <- rmultinom(1, size = new_totals[i],
                              prob = counts[,i] / sum(counts[,i]))
    }
  }

  counts <- data.matrix(counts)
  colnames(counts) <- NULL
  rownames(counts) <- NULL
  
  parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
  return(parsed_obj)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#
#   Hashimshony et al. (2016)
#
#   No longer in use!
#   Very little detectable differential abundance in "absolute" counts here
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Parse the Hashimshony et al. (2016) mouse fibroblast data set
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @import edgeR
#' @export
parse_Hashimshony <- function(absolute = TRUE, use_cpm = FALSE) {
  # Note that 0 = GFP negative cells (quiescent)
  #           1 = GFP positive cells (cycling)
  #           0.3 = GFP "weak"
  #           0.5 = GFP "mixed"
  file_dir <- file.path("data", "Hashimshony_2016")
  
  # Parse data and assignments
  data <- read.table(file.path(file_dir,
                               "GSE78779_Expression_C1_96_cells.txt"),
                     sep = "\t",
                     header = TRUE)
  
  groups <- sapply(colnames(data)[2:ncol(data)], function(x) {
    strsplit(x, "_")[[1]][2]
  })
  # Eliminate samples with missing/ambiguous cell cycle markers
  retain_idx <- which(!(groups %in% c(".", "n.a.")))
  data <- data[,c(1, retain_idx+1)]
  groups <- groups[retain_idx]
  groups <- factor(unname(groups), levels = c("0", "0.3", "0.5", "1"))
  
  spike_idx <- which(sapply(data$Sample, function(x) str_detect(x, "^ERCC-")))
  data <- data[,2:ncol(data)]

  if(use_cpm) {
    data <- apply(data, 2, function(x) x/sum(x))
    data <- round(data*1e06)
  }
  
  # Here total abundances look pretty informative. Use:
  #   gapdh_idx <- which(data$Sample == "ENSMUSG00000057666")
  # which returns 16179 to compare Gadph abundance to totals
  sf <- compute_sf(data[spike_idx,])
  counts <- data[-spike_idx,]

  if(absolute) {
    for(i in 1:ncol(counts)) {
      counts[,i] <- counts[,i] / sf[i]
    }
  } else {
    set.seed(101)
    new_totals <- sample(round(colSums(counts)))
    for(i in 1:ncol(counts)) {
      counts[,i] <- rmultinom(1, size = new_totals[i],
                              prob = counts[,i] / sum(counts[,i]))
    }
  }
  counts <- data.matrix(counts)
  colnames(counts) <- NULL
  rownames(counts) <- NULL
  
  # Subset to the SC samples
  A_idx <- which(groups == "0")
  B_idx <- which(groups == "1")
  counts <- counts[,c(A_idx, B_idx)]
  groups <- c(rep("0", length(A_idx)), rep("1", length(B_idx)))
  
  parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
  return(parsed_obj)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#
#   Kimmerling et al. (2018)
#
#   No longer in use!
#   Very little detectable differential abundance in "absolute" counts here
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Parse the Kimmerling et al. (2018) fibroblast data set
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @export
parse_Kimmerling <- function(absolute = TRUE, use_spike_ins = FALSE, use_cpm = FALSE) {
  file_dir <- file.path("data", "Kimmerling_2018")
  
  data <- read.table(file.path(file_dir, "fl5_serial_rsem3.txt"),
                     sep = "\t",
                     header = TRUE,
                     row.names = 1)
  meta <- read.table(file.path(file_dir, "qc_fl5_serial3.txt"),
                     sep = "\t",
                     header = TRUE,
                     row.names = 1)
  
  # Column names in `data` appear to match the rownames in `meta`
  # Rows (samples) in `meta` are alphabetically ordered; arrange `data` columns
  # to match
  reorder <- order(colnames(data))
  data <- data[,reorder]
  
  # QC as in the author's analysis
  qc_samples <- meta$mass >= 0 & meta$gene_count >= 500
  # qc_samples <- meta$mass >= 0
  
  data <- data[,qc_samples]
  meta <- meta[qc_samples,]
  
  # The rationale for this squick and dirty renormalization is that I want a
  # size factor generated from the distribution of mass, centered at 1, with
  # a SD of around 0.1
  # sf <- scale(meta$mass/(sd(meta$mass)*10), scale = F) + 1
  sf <- compute_sf(meta$mass)
  counts <- data
  
  if(use_cpm) {
    counts <- apply(counts, 2, function(x) x/sum(x))
    counts <- round(counts*1e06)
  }
  
  if(absolute) {
    for(i in 1:ncol(counts)) {
      # Multiply; library size should have a positive association with mass!
      counts[,i] <- counts[,i] * sf[i]
    }
  } else {
    set.seed(101)
    new_totals <- sample(round(colSums(counts)))
    for(i in 1:ncol(counts)) {
      counts[,i] <- rmultinom(1, size = new_totals[i],
                              prob = counts[,i] / sum(counts[,i]))
    }
  }
  
  # We'll separate groups into highest and lowest thirds by mass
  thresholds <- quantile(meta$mass, probs = c(1/3, 2/3))
  groups <- numeric(ncol(counts))
  groups[meta$mass < thresholds[1]] <- 1
  groups[meta$mass > thresholds[2]] <- 2
  
  counts <- counts[,groups != 0]
  meta <- meta[groups != 0,]
  groups <- groups[groups != 0]
  
  # Reorder the data and metadata
  reorder <- order(groups)
  counts <- counts[,reorder]
  meta <- meta[reorder,]
  groups <- groups[reorder]
  
  counts <- data.matrix(counts)
  colnames(counts) <- NULL
  rownames(counts) <- NULL
  
  groups <- factor(groups, levels = c("1", "2"))
  levels(groups) <- c("low_mass", "high_mass")
  
  parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
  return(parsed_obj)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#
#   TCGA ESCA RNA-seq data
#
#   No longer in use!
#   Very little fold change between conditions here
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Parse the TCGA ESCA RNA-seq data
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @export
parse_ESCA <- function(absolute = TRUE) {
  file_dir <- file.path("data", "TCGA_ESCA")

  # Pull SCC labels from Campbell et al. supplement
  sccs <- read.delim(file.path(file_dir,
                               "cell_reports_squamous_TCGA_paper_sup1_mmc2_table_S1M.txt"))
  
  # Pull ESCA log2 RPKM RNA-seq data
  esca <- read.delim(file.path(file_dir,
                               "Human__TCGA_ESCA__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct"))
  colnames(esca)[1] <- "gene_name"
  
  spikein_idx <- unname(which(sapply(esca$gene_name, function(x) {
    x == "GAPDH"
  })))

  # Convert log2 pseudocounts back to counts  
  esca <- 2**esca[,2:ncol(esca)] - 1

  # Find spike-ins and separate into spike-in and all other sequence data sets
  esca_spike <- esca[spikein_idx,]
  esca <- esca[-spikein_idx,]
  
  # Label samples with tumor subtype
  barcodes <- unname(sapply(colnames(esca), function(x) {
    str_replace_all(x, "\\.", "-")
  }))
  scc_flag <- rep("Other", ncol(esca))
  scc_flag[barcodes %in% sccs$Patient.Barcodes] <- "SCC"
  
  # Plot omitted here but there's an overall higher spike-in abundance in the
  # SCC samples
  if(absolute) {
    sf <- compute_sf(esca_spike)
    esca <- esca / sf
  } else {
    set.seed(101)
    new_totals <- sample(round(colSums(esca)))
    for(i in 1:ncol(esca)) {
      esca[,i] <- rmultinom(1, size = new_totals[i],
                              prob = esca[,i] / sum(esca[,i]))
    }
  }
  groups <- scc_flag
  tax = NULL
  rownames(esca) <- NULL
  colnames(esca) <- NULL
  esca <- data.matrix(esca)
  
  esca <- cbind(esca[,groups == "Other"], esca[,groups == "SCC"])
  groups <- sort(groups)
  
  parsed_obj <- list(counts = esca, groups = groups, tax = NULL)
  return(parsed_obj)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#
#   Gruen et al. (2014)
#
#   No longer in use!
#   Very little fold change between conditions here
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Parse the Gruen et al. (2014) single-cell cell line data set
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @export
parse_Gruen <- function(absolute = TRUE, use_cpm = FALSE) {
  counts <- read.table(file.path("data",
                                 "Gruen_2014",
                                 "GSE54695_data_transcript_counts.txt"),
                       header = TRUE)
  gene_IDs <- counts$GENENAME
  counts <- counts[,2:ncol(counts)]
  
  # First remove samples with low counts
  counts <- counts[,colSums(counts) > 5000]
  
  platform <- unname(sapply(colnames(counts), function(x) {
    str_split(x, "_")[[1]][1]
  }))
  condition <- unname(sapply(colnames(counts), function(x) {
    str_split(x, "_")[[1]][2]
  }))
  
  # Pull spike-in sequences
  ref_idx <- which(sapply(gene_IDs, function(x) str_detect(x, "^ERCC-")))
  spike_counts <- counts[ref_idx,]
  sf <- unname(compute_sf(spike_counts))
  
  if(use_cpm) {
    counts <- apply(counts, 2, function(x) x/sum(x))
    counts <- round(counts*1e06)
  }
  
  counts <- counts[-ref_idx,]
  
  if(absolute) {
    for(j in 1:ncol(counts)) {
      counts[,j] <- counts[,j]/sf[j]
    }
  } else {
    # Shuffle observed abundances
    set.seed(1000)
    new_totals <- sample(round(colSums(counts)))
    for(i in 1:ncol(counts)) {
      counts[,i] <- rmultinom(1, size = new_totals[i],
                              prob = counts[,i] / sum(counts[,i]))
    }
  }
  rownames(counts) <- NULL
  colnames(counts) <- NULL
  counts <- data.matrix(counts)
  
  # Subset to the SC samples
  A_idx <- which(platform == "SC" & condition == "2i")
  B_idx <- which(platform == "SC" & condition == "serum")
  counts <- counts[,c(A_idx, B_idx)]
  groups <- c(rep("A", length(A_idx)), rep("B", length(B_idx)))
  
  parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
  return(parsed_obj)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#
#   Ferreira et al. (2014)
#
#   No longer in use!
#   Very little detectable differential abundance in "absolute" counts here
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Parse the Ferreira et al. (2014) zebrafish bulk RNA-seq data set
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @export
parse_Ferreira <- function(absolute = TRUE) {
  counts_A1 <- read.table(file.path("data",
                                    "Ferreira_2014",
                                    "GSM1289430_Control_1_rep1.txt"))
  counts_A2 <- read.table(file.path("data",
                                    "Ferreira_2014",
                                    "GSM1289431_Control_2_rep1.txt"))
  counts_A3 <- read.table(file.path("data",
                                    "Ferreira_2014",
                                    "GSM1289432_Control_3_rep1.txt"))
  counts_B1 <- read.table(file.path("data",
                                    "Ferreira_2014",
                                    "GSM1289433_Treated_1_rep1.txt"))
  counts_B2 <- read.table(file.path("data",
                                    "Ferreira_2014",
                                    "GSM1289434_Treated_2_rep1.txt"))
  counts_B3 <- read.table(file.path("data",
                                    "Ferreira_2014",
                                    "GSM1289435_Treated_3_rep1.txt"))
  
  clip_row <- which(counts_A1$V1 == "no_feature")
  counts_A1 <- counts_A1[1:(clip_row-1),]
  counts_A2 <- counts_A2[1:(clip_row-1),]
  counts_A3 <- counts_A3[1:(clip_row-1),]
  counts_B1 <- counts_B1[1:(clip_row-1),]
  counts_B2 <- counts_B2[1:(clip_row-1),]
  counts_B3 <- counts_B3[1:(clip_row-1),]
  
  gene_IDs <- counts_A1$V1
  counts_A <- cbind(counts_A1[,2], counts_A2[,2], counts_A3[,2])
  counts_B <- cbind(counts_B1[,2], counts_B2[,2], counts_B3[,2])
  
  counts <- cbind(counts_A, counts_B)
  spike_idx <- which(sapply(gene_IDs, function(x) str_detect(x, "^ERCC-")))
  spike_counts <- counts[spike_idx,]
  sf <- compute_sf(spike_counts)
  counts <- counts[-spike_idx,]
  
  if(absolute) {
    for(j in 1:ncol(counts)) {
      counts[,j] <- counts[,j]/sf[j]
    }
  } else {
    # Shuffle observed abundances
    set.seed(101)
    new_totals <- sample(round(colSums(counts)))
    plot(new_totals)
    for(i in 1:ncol(counts)) {
      counts[,i] <- rmultinom(1, size = new_totals[i],
                              prob = counts[,i] / sum(counts[,i]))
    }
  }
  
  groups <- c(rep("control", 3), rep("treated", 3))
  parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
  return(parsed_obj)
}
