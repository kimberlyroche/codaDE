#' Parse the Barlow et al. (2020) 16S data
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @export
parse_Barlow <- function(absolute = TRUE) {
  if(absolute) {
    parsed_data_fn <- file.path("data", "absolute_Barlow.rds")
  } else {
    parsed_data_fn <- file.path("data", "relative_Barlow.rds")
  }
  if(file.exists(parsed_data_fn)) {
    return(readRDS(parsed_data_fn))
  } else {
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
    counts <- as.matrix(counts)
    colnames(counts) <- NULL
    rownames(counts) <- NULL
    
    if(!absolute) {
      # Note relative_counts are scaled to 100; scale to CPM
      counts <- counts * (1e6 / sum(counts[1,]))
    }
    
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
    
    groups <- factor(labels, levels = c("control", "keto"))
    
    # The scaling does weird things to this data set. The minimum *observed* 
    # (non-zero) abundance is huge giving a gap between observed and unobserved 
    # features that is enormous. I'm scaling down all the counts by this minimum 
    # observed abundance for now.
    min_observed <- min(counts[which(counts > 0, arr.ind = TRUE)])
    counts <- counts / min_observed
    counts <- t(counts)
    
    parsed_obj <- list(counts = counts, groups = groups, tax = tax)
    saveRDS(parsed_obj, file = parsed_data_fn)
    return(parsed_obj)
  }
}

#' Parse the Morton et al. (2019) 16S data
#'
#' @return named list of counts and group labels
#' @export
parse_Morton <- function() {
  parsed_data_fn <- file.path("data", "absolute_Morton.rds")
  if(file.exists(parsed_data_fn)) {
    return(readRDS(parsed_data_fn))
  } else {
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
    
    parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
    saveRDS(parsed_obj, file = parsed_data_fn)
    return(parsed_obj)
  }
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

#' Parse the Athanasiadou et al. (2021) bulk RNA-seq data
#'
#' @param which_data options are "ciona" or "yeast"
#' @return named list of counts and group labels
#' @export
parse_Athanasiadou <- function(which_data = "ciona") {
  which_data <- tolower(which_data)
  if(!(which_data %in% c("ciona", "yeast"))) {
    stop("Invalid data type!")
  }
  parsed_data_fn <- file.path("data",
                              paste0("absolute_Athanasiadou_",
                                     which_data,
                                     ".rds"))
  if(file.exists(parsed_data_fn)) {
    return(readRDS(parsed_data_fn))
  } else {
    if(which_data == "yeast") {
      file_dir <- file.path("data", "Athanasiadou_2021", "S1CodeandData")
      
      groups <- as.factor(rep(c("C12", "C20", "C30"), rep(3,3)))
      counts <- normalize_Athanasiadou(file_dir = file_dir,
                                       SI_file = "YmatSIyeast.txt",
                                       RNA_file = "YmatRNAyeast.txt",
                                       ERCC_annot = "ERCCyeast.txt",
                                       groups = groups)
      
      parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
      saveRDS(parsed_obj, file = parsed_data_fn)
      return(parsed_obj)
    } else {
      # Ciona robusta
      file_dir <- file.path("data", "Athanasiadou_2021", "S1CodeandData")
      
      groups <- as.factor(rep(c("lacz", "dnfgfr", "camras"), rep(3,3)))
      counts <- normalize_Athanasiadou(file_dir = file_dir,
                                       SI_file = "YmatSIciona.txt",
                                       RNA_file = "YmatRNAciona.txt",
                                       ERCC_annot = "ERCCciona.txt",
                                       groups = groups)
      
      baseline_counts <- counts[,which(groups == "lacz")]
      absent_tax_idx <- which(rowSums(baseline_counts) == 0)
      baseline_counts <- baseline_counts[,-absent_tax_idx]
      mean_ciona <- mean(baseline_counts)
      
      # Ciona counts are super low relative to (abundant) spike-in
      # Scale these up to the mean counts in the S. cerevisieae experiment
      mean_yeast <- 7.921205 # I've pre-calculated this, as above
      scale <- mean_yeast / mean_ciona
      counts <- counts * scale
      
      parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
      saveRDS(parsed_obj, file = parsed_data_fn)
      return(parsed_obj)
    }
  }
}
