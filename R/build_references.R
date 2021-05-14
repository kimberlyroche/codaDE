#' Generate simulated differential expression for two conditions
#'
#' @param p number of features (genes, taxa) to simulate
#' @param log_mean log mean for condition 1
#' @param log_var log variance for condition 1
#' @param log_noise_var log variance for noise component added to condition 1
#' to give condition 2
#' @param base_correlation this is the base correlation matrix for simulated
#' features in log space
#' @param concentration concentration parameter for sampling the correlation 
#' across feature; if base_correlation is NULL and concentration = LARGE, features
#' are effectively independent
#' @param save_name optional tag to append to saved data file
#' @return NULL
#' @import matrixsampling
#' @import MASS
#' @export
build_simulated_reference <- function(p = 1000, log_mean = 0, log_var = 2,
                                      log_noise_var = 1, base_correlation = NULL,
                                      concentration = 1e6, save_name = NULL) {
  if(is.null(base_correlation)) {
    base_correlation <- diag(p)
  }
  if(concentration < p + 2) {
    stop("Invalid concentration specified!")
  }
  # log_counts1 <- rnorm(p, log_mean, log_var)
  # log_counts2 <- log_counts1 + rnorm(p, 0, log_noise_var)
  K <- cov2cor(rinvwishart(1, concentration, base_correlation)[,,1])
  
  log_counts1 <- mvrnorm(1, rep(log_mean, p), diag(p)*(log_var^2))
  log_perturbation <- mvrnorm(1, rep(0, p), K*(log_noise_var^2))
  log_counts2 <- log_counts1 + log_perturbation

  counts1 <- exp(log_counts1)
  counts2 <- exp(log_counts2)
  
  sim_obj <- list(log_cond1 = log_counts1, log_cond2 = log_counts2,
                  cond1 = counts1, cond2 = counts2,
                  log_perturbation = log_perturbation, correlation_matrix = K)
  
  if(is.null(save_name)) {
    return(sim_obj)
  } else {
    # save_file <- ifelse(is.null(save_name),
    #                     "DE_reference_simulated.rds",
    #                     paste0("DE_reference_",save_name,".rds"))
    save_file <- paste0("DE_reference_",save_name,".rds")
    saveRDS(sim_obj,
            file = file.path("data", save_file))
  }
}

#' Generate differential expression reference for Barlow et al. (2020) 16S data
#'
#' @return NULL
#' @export
build_Barlow_reference <- function() {
  file_dir <- file.path("data", "Barlow_2020")
  
  data <- read.table(file.path(file_dir, "Absolute_Abundance_Table.csv"),
                     header = TRUE, sep = ",", stringsAsFactors = FALSE)

  # Clean up data
  md_labels <- c("Diet", "Site", "Day", "mouse", "Cage")
  metadata <- data[,md_labels]
  counts <- data[,!(colnames(data) %in% md_labels)]
  counts <- counts[,2:ncol(counts)]
  tax <- colnames(counts)
  counts <- as.matrix(counts)
  colnames(counts) <- NULL
  rownames(counts) <- NULL

  # `counts` is initially 103 samples x 142 taxa

  # Pull stool samples from days 4, 7, 10
  keto_idx <- which(metadata$Diet == "Keto" & metadata$Day %in% c(4,7,10) & metadata$Site == "Stool")
  ctrl_idx <- which(metadata$Diet == "Control" & metadata$Day %in% c(4,7,10) & metadata$Site == "Stool")

  days <- metadata[c(ctrl_idx,keto_idx),]$Day
  
  counts <- rbind(counts[ctrl_idx,], counts[keto_idx,])
  labels <- c(rep("control", length(ctrl_idx)), rep("keto", length(keto_idx)))

  # Eliminate all-zero taxa
  absent_tax_idx <- which(colSums(counts) == 0)
  counts <- counts[,-absent_tax_idx]
  tax <- tax[-absent_tax_idx]

  groups = factor(labels, levels = c("control", "keto"))

  # The scaling does weird things to this data set. The minimum *observed* 
  # (non-zero) abundance is huge giving a gap between observed and unobserved 
  # features that is enormous. I'm scaling down all the counts by this minimum 
  # observed abundance for now.
  min_observed <- min(counts[which(counts > 0, arr.ind = TRUE)])
  counts <- counts / min_observed
  counts <- t(counts)
  
  saveRDS(list(counts = counts, groups = groups, tax = tax),
         file = file.path("data", "absolute_Barlow.rds"))
  
  # Build DE model
  counts1 <- counts[,groups == "control" & days == 10]
  counts2 <- counts[,groups == "keto" & days == 10]
  
  # Arbitrarily "match" a subject on day 10 in condition 1 with a subject on
  # day 10 in condition 2
  
  counts1 <- c(counts1)
  counts2 <- c(counts2)
  
  saveRDS(list(cond1 = counts1, cond2 = counts2),
          file = file.path("data", "DE_reference_Barlow.rds"))
  
  return(NULL)
}

#' Generate differential expression reference for Morton et al. (2019) 16S data
#'
#' @return NULL
#' @export
build_Morton_reference <- function() {
  # TBD - move bipartite graph code
  
  file_dir <- file.path("data", "Morton_2019", "Github")

  otu_table <- read.delim(file.path(file_dir, "oral_trimmed_deblur.txt"),
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  metadata <- read.delim(file.path(file_dir, "oral_trimmed_metadata.txt"),
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  # "flowcount" is calculated this way in the paper
  # (See their Python notebook for Fig. 2)
  flowcount <- (metadata$flow.cell.5min.1 + metadata$flow.cell.5min.2) / 2
  event_labels <- metadata$brushing_event
  subject_labels <- metadata$HostSubject
  timepoint_labels <- metadata$`Timepoint.`
  timepoint_labels[timepoint_labels %in% c(1,2)] <- "_AM"
  timepoint_labels[timepoint_labels %in% c(7,8,9)] <- "_PM"
  sample_labels <- metadata$SampleID

  # Compute relative abundances
  tax_legend <- otu_table$OTU_ID
  props <- otu_table[,2:ncol(otu_table)]
  props <- apply(props, 2, function(x) x/sum(x))

  # Scale these up to absolute abundances
  # `props` is taxa x samples
  counts <- props
  for(i in 1:ncol(counts)) {
    counts[,i] <- counts[,i]*flowcount[i]
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

  # Scale by minimum observed abundance as I did with Barlow et al.
  min_observed <- min(counts[which(counts > 0, arr.ind = TRUE)])
  counts <- counts / min_observed
  
  saveRDS(list(counts = counts, groups = groups, tax = NULL),
          file = file.path("data", "absolute_Morton.rds"))

  # Build DE model
  counts1 <- c(counts[,1:length(before_events)])
  counts2 <- c(counts[,(length(before_events)+1):ncol(counts)])

  saveRDS(list(cond1 = counts1, cond2 = counts2),
          file = file.path("data", "DE_reference_Morton.rds"))
  
  return(NULL)
}

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

#' Generate differential expression reference for Athanasiadou et al. (2021) 
#' bulk RNA-seq data
#'
#' @return named list of observed counts in conditions 1 and 2
#' @export
build_Athanasiadou_reference <- function() {
  # TBD - allow choice of baseline v. treatment
  file_dir <- file.path("data", "Athanasiadou_2021", "S1CodeandData")
  
  groups <- as.factor(rep(c("C12", "C20", "C30"), rep(3,3)))
  counts <- normalize_Athanasiadou(SI_file = "YmatSIyeast.txt",
                                   RNA_file = "YmatRNAyeast.txt",
                                   ERCC_annot = "ERCCyeast.txt",
                                   groups = groups)
  
  # Histograms!
  baseline_counts <- counts[,which(groups == "C12")]
  absent_tax_idx <- which(rowSums(baseline_counts) == 0)
  baseline_counts <- baseline_counts[,-absent_tax_idx]
  mean_yeast <- mean(baseline_counts)
  
  saveRDS(list(counts = counts, groups = groups),
          file = file.path("data", "absolute_Athanasiadou_yeast.rds"))
  
  # Build DE model (arbitrarily matching replicates)
  counts1 <- c(counts[,groups == "C12"])
  counts2 <- c(counts[,groups == "C30"])
  
  saveRDS(list(cond1 = counts1, cond2 = counts2),
          file = file.path("data", "DE_reference_Athanasiadou_yeast.rds"))
  
  groups <- as.factor(rep(c("lacz", "dnfgfr", "camras"), rep(3,3)))
  counts <- normalize_Athanasiadou(SI_file = "YmatSIciona.txt",
                                   RNA_file = "YmatRNAciona.txt",
                                   ERCC_annot = "ERCCciona.txt",
                                   groups = groups)
  
  baseline_counts <- counts[,which(groups == "lacz")]
  absent_tax_idx <- which(rowSums(baseline_counts) == 0)
  baseline_counts <- baseline_counts[,-absent_tax_idx]
  mean_ciona <- mean(baseline_counts)

  # Ciona counts are super low relative to (abundant) spike-in
  # Scale these up to the mean counts in the S. cerevisieae experiment
  scale <- mean_yeast / mean_ciona
  counts <- counts * scale
  
  saveRDS(list(counts = counts, groups = groups),
          file = file.path("data", "absolute_Athanasiadou_ciona.rds"))
  
  # Build DE model
  counts1 <- c(counts[,groups == "lacz"])
  counts2 <- c(counts[,groups == "dnfgfr"])
  
  saveRDS(list(cond1 = counts1, cond2 = counts2),
          file = file.path("data", "DE_reference_Athanasiadou_ciona.rds"))
}
