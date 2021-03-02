#' Generate simulated differential expression for two conditions
#'
#' @param p number of features (genes, taxa) to simulate
#' @param log_var log variance for condition 1
#' @return log_noise_var log variance for noise component added to condition 1
#' to give condition 2
#' @export
build_simulated_reference <- function(p = 1000, log_var = 2, log_noise_var = 1) {
  log_counts1 <- rnorm(p, 0, log_var)
  log_counts2 <- log_counts1 + rnorm(p, 0, log_noise_var)
  counts1 <- exp(log_counts1)
  counts2 <- exp(log_counts2)
  saveRDS(list(cond1 = counts1, cond2 = counts2),
          file = file.path("data", "DE_reference_simulated.rds"))
}

#' Generate differential expression reference for Barlow et al. (2020) 16S data
#'
#' @return named list of observed counts in conditions 1 and 2
#' @export
build_Barlow_reference <- function() {
  file_dir <- file.path("data", "Barlow_2020")
  
  data <- read.table(file.path(file_dir, "Absolute_Abundance_Table.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)

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

  # # Histograms!
  # plot_counts <- counts[which(groups == "control"),]
  # absent_tax_idx <- which(colSums(plot_counts) == 0)
  # plot_counts <- plot_counts[,-absent_tax_idx]
  # 
  # ggplot(data.frame(x = log(c(plot_counts) + 1)), aes(x)) +
  #   geom_histogram(color = "white")
  # ggsave("hist_Barlow_orig.png", units = "in", dpi = 100, height = 5, width = 6)
  # 
  # ggplot(data.frame(x = log(c(plot_counts)/min(plot_counts[which(plot_counts > 0, arr.ind = TRUE)]) + 1)), aes(x)) +
  #   geom_histogram(color = "white")
  # ggsave("hist_Barlow_scaled.png", units = "in", dpi = 100, height = 5, width = 6)
  # 
  # The scaling does weird things to this data set. The minimum *observed* (non-zero) abundance is huge
  # giving a gap between observed and unobserved features that is enormous.
  # I'm scaling down all the counts by this minimum observed abundance for now. (2/16/2021)
  min_observed <- min(counts[which(counts > 0, arr.ind = TRUE)])
  counts <- counts / min_observed
  counts <- t(counts)
  
  saveRDS(list(counts = counts, groups = groups, tax = tax),
         file = file.path(file_dir, "absolute_Barlow.rds"))
  
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
#' @return named list of observed counts in conditions 1 and 2
#' @export
build_Morton_reference <- function() {
  # TBD - move bipartite graph code
  
  file_dir <- file.path("data", "Morton_2019", "Github")

  otu_table <- read.delim(file.path(file_dir, "oral_trimmed_deblur.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  metadata <- read.delim(file.path(file_dir, "oral_trimmed_metadata.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  # "flowcount" is calculated this way in the paper (see Python notebook for Fig. 2)
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

  # Re-order samples by condition (event before/after) and subject (where most subjects have a morning and evening brushing event)
  before_events <- which(event_labels == "before")
  before_subjects <- paste0(subject_labels[before_events], timepoint_labels[before_events])
  after_events <- which(event_labels == "after")
  after_subjects <- paste0(subject_labels[after_events], timepoint_labels[after_events])
  counts <- counts[,c(before_events, after_events)]
  subjects <- factor(c(before_subjects, after_subjects))

  groups <- factor(c(rep("before", length(before_events)), rep("after", length(after_events))),
                   levels = c("before", "after"))

  # Scale by minimum observed abundance as I did with Barlow et al.
  min_observed <- min(counts[which(counts > 0, arr.ind = TRUE)])
  counts <- counts / min_observed
  
  # # Bipartite graph
  # # features <- sample(1:ncol(counts), size = 5)
  # plot_data <- data.frame(before_log_count = c(), after_log_count = c(), rank = c())
  # before_log_counts <- log(unname(counts[,groups == "before"]) + 1)
  # qq <- c(sum(before_log_counts == 0)/(nrow(before_log_counts)*ncol(before_log_counts)))
  # qq <- quantile(before_log_counts, probs = seq(from = qq, to = 1, length.out = 5))
  # breaks <- c(-Inf, qq, Inf)
  # sampled_subjects <- sample(before_subjects, size = 5)
  # for(feat in 1:nrow(counts)) {
  #   for(subj in sampled_subjects) {
  #     before_val <- log(unname(counts[feat,subjects == subj & groups == "before"]) + 1)
  #     after_val <- log(unname(counts[feat,subjects == subj & groups == "after"]) + 1)
  #     rank <- as.numeric(cut(before_val, breaks = breaks))
  #     plot_data <- rbind(plot_data, data.frame(before_log_count = before_val,
  #                                              after_log_count = after_val,
  #                                              rank = rank))
  #   }
  # }
  # plot_data$rank <- as.factor(plot_data$rank)
  # 
  # ggplot(plot_data, aes(x = 0, xend = 1, y = before_log_count, yend = after_log_count, color = rank)) +
  #   geom_segment(size = 1.2, alpha = 0.3)
  
  saveRDS(list(counts = counts, groups = groups, tax = NULL),
          file = file.path("data", "absolute_Morton.rds"))

  # Build DE model
  counts1 <- c(counts[,1:length(before_events)])
  counts2 <- c(counts[,(length(before_events)+1):ncol(counts)])

  saveRDS(list(cond1 = counts1, cond2 = counts2),
          file = file.path("data", "DE_reference_Morton.rds"))
  
  return(NULL)
}

# Utility function -- applies the ERCC normalization code from Athanasiadou et al. (2021)
normalize_Athanasiadou <- function(SI_file, RNA_file, ERCC_annot, groups) {
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
  
  # From RMarkdown comments: "Compute spike-in proportions. Choose as the reference spike-in,
  #   the one with the largest proportion; record its proportion and attomoles; and compute the
  # library-specific $\nu_j$ calibration factors"
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

#' Generate differential expression reference for Athanasiadou et al. (2021) bulk RNA-seq data
#'
#' @return named list of observed counts in conditions 1 and 2
#' @export
build_Athanasiadou_reference <- function() {
  # TBD - allow choice of baseline v. treatment
  #       move bipartite graph code
  
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
  
  # Bipartite graph
  # label1 <- "lacz"
  # label1 <- "camras"
  # label2 <- "dnfgfr"
  # plot_data <- data.frame(before_log_count = c(), after_log_count = c(), rank = c())
  # before_log_counts <- log(unname(counts[,groups == label1]) + 1)
  # # qq <- c(sum(before_log_counts == 0)/(nrow(before_log_counts)*ncol(before_log_counts)))
  # # qq <- quantile(before_log_counts, probs = seq(from = qq, to = 1, length.out = 5))
  # qq <- quantile(before_log_counts, probs = seq(from = 0, to = 1, length.out = 6))
  # breaks <- c(-Inf, qq[2:(length(qq)-1)], Inf)
  # sampled_features <- sample(1:nrow(counts), size = 200)
  # for(feat in sampled_features) {
  #   for(samp in 1:3) {
  #     before_val <- log(unname(counts[feat,groups == label1][samp]) + 1)
  #     after_val <- log(unname(counts[feat,groups == label2][samp]) + 1)
  #     rank <- as.numeric(cut(before_val, breaks = breaks))
  #     plot_data <- rbind(plot_data, data.frame(before_log_count = before_val,
  #                                              after_log_count = after_val,
  #                                              rank = rank))
  #   }
  # }
  # plot_data$rank <- as.factor(plot_data$rank)
  # 
  # ggplot(plot_data, aes(x = 0, xend = 1, y = before_log_count, yend = after_log_count, color = rank)) +
  #   geom_segment(size = 1.2, alpha = 0.3)
  
  saveRDS(list(counts = counts, groups = groups),
          file = file.path("data", "absolute_Athanasiadou_ciona.rds"))
  
  # Build DE model
  counts1 <- c(counts[,groups == "lacz"])
  counts2 <- c(counts[,groups == "dnfgfr"])
  
  saveRDS(list(cond1 = counts1, cond2 = counts2),
          file = file.path("data", "DE_reference_Athanasiadou_ciona.rds"))
}
