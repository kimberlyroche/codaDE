#' Compute reference derived per-sample scale factor
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

#' Parse the Barlow et al. (2020) 16S data
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @export
parse_Barlow <- function(absolute = TRUE) {
  # if(absolute) {
  #   parsed_data_fn <- file.path("data", "absolute_Barlow.rds")
  # } else {
  #   parsed_data_fn <- file.path("data", "relative_Barlow.rds")
  # }
  # if(file.exists(parsed_data_fn)) {
  #   return(readRDS(parsed_data_fn))
  # } else {
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
    # saveRDS(parsed_obj, file = parsed_data_fn)
    return(parsed_obj)
  # }
}

#' Parse the Morton et al. (2019) 16S data
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @export
parse_Morton <- function(absolute = TRUE) {
  # parsed_data_fn <- file.path("data", "absolute_Morton.rds")
  # if(file.exists(parsed_data_fn)) {
  #   return(readRDS(parsed_data_fn))
  # } else {
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
    
    parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
    # saveRDS(parsed_obj, file = parsed_data_fn)
    return(parsed_obj)
  # }
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
  # parsed_data_fn <- file.path("data",
  #                             paste0("absolute_Athanasiadou_",
  #                                    which_data,
  #                                    ".rds"))
  # if(file.exists(parsed_data_fn)) {
  #   return(readRDS(parsed_data_fn))
  # } else {
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
      
      parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
      # saveRDS(parsed_obj, file = parsed_data_fn)
      return(parsed_obj)
    }
  # }
}

#' Parse the Song et al. (2021) nCounter data
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @import edgeR
#' @export
parse_Song <- function(absolute = TRUE) {
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
  
  # Use edgeR (for now) to compute the size factor
  # sf <- calcNormFactors(mrna[ref_idx,]) # columns assumed to be samples
  counts <- as.matrix(mrna[-ref_idx,])
  if(!absolute) {
    # Shuffle observed abundances
    set.seed(1000)
    new_totals <- sample(round(colSums(counts)))
    for(i in 1:ncol(counts)) {
      counts[,i] <- rmultinom(1, size = new_totals[i], prob = counts[,i] / sum(counts[,i]))
    }
  }
  # Spike-in renormalization
  # if(absolute) {
  #   for(j in 1:ncol(counts)) {
  #     counts[,j] <- counts[,j] / sf[j]
  #   }
  # }
  colnames(counts) <- NULL
  rownames(counts) <- NULL
  
  parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
  return(parsed_obj)
}

#' Parse the Muraro et al. (2016) single-cell pancreas data set
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @import edgeR
#' @export
parse_Muraro <- function(absolute = TRUE) {
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
  c1 <- data.frame(expr = unlist(unname(data[gcg_idx,])),
                   cluster = mapping$cluster[!is.na(mapping$cluster)]) %>%
    group_by(cluster) %>%
    summarize(mean_expr = mean(expr)) %>%
    arrange(desc(mean_expr)) %>%
    top_n(1) %>%
    pull(cluster)
  c2 <- data.frame(expr = unlist(unname(data[ins_idx,])),
                   cluster = mapping$cluster[!is.na(mapping$cluster)]) %>%
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
  
  # Pull spike-in sequences
  spikein_seqs <- which(sapply(rownames(data), function(x) str_detect(x, "^ERCC-\\d+")))
  spikein_counts <- cbind(data_orig[spikein_seqs,
                                    mapping %>% filter(!is.na(cluster)) %>% filter(cluster %in% c(c1)) %>% pull(idx)],
                          data_orig[spikein_seqs,
                                    mapping %>% filter(!is.na(cluster)) %>% filter(cluster %in% c(c2)) %>% pull(idx)])
  
  # Use edgeR (for now) to compute the size factor
  # sf <- calcNormFactors(spikein_counts) # columns assumed to be samples
  sf <- compute_sf(spikein_counts)
  
  counts <- counts[-spikein_seqs,]
  if(absolute) {
    for(j in 1:ncol(counts)) {
      counts[,j] <- counts[,j] / sf[j]
    }
  }
  counts <- as.matrix(counts)
  colnames(counts) <- NULL
  rownames(counts) <- NULL
  
  parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
  return(parsed_obj)
}

#' Parse the Monaco et al. (2019) immune cell data set
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @import edgeR
#' @export
parse_Monaco <- function(absolute = TRUE) {
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

  # Plot spike-in relative abundance by type
  groups <- sapply(colnames(data), function(x) {
    pieces <- strsplit(x, "_")[[1]]
    paste0(pieces[2:length(pieces)], collapse = "_")
  })
  
  # Compute size factor
  # sf <- unname(unlist(calcNormFactors(data[spike_idx,]))) # columns assumed to be samples
  sf <- compute_sf(data[spike_idx,])
  
  # Re-normalize
  counts <- data
  if(absolute) {
    for(i in 1:ncol(counts)) {
      counts[,i] <- counts[,i] / sf[i]
    }
  }
  counts <- as.matrix(counts)
  colnames(counts) <- NULL
  rownames(counts) <- NULL
  
  plot(colSums(counts))
  
  parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
  return(parsed_obj)
}

#' Parse the Hashimshony et al. (2016) mouse fibroblast data set
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @import edgeR
#' @export
parse_Hashimshony <- function(absolute = TRUE) {
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
  
  # Here total abundances look pretty informative. Use:
  #   gapdh_idx <- which(data$Sample == "ENSMUSG00000057666")
  # which returns 16179
  
  # Re-normalize
  counts <- data
  if(absolute) {
    # Compute size factor from spike-ins
    # sf <- unname(unlist(calcNormFactors(counts[spike_idx,]))) # columns assumed to be samples
    sf <- compute_sf(counts[spike_idx,])
    for(i in 1:ncol(counts)) {
      counts[,i] <- counts[,i] / sf[i]
    }
  }
  counts <- as.matrix(counts)
  colnames(counts) <- NULL
  rownames(counts) <- NULL
  
  parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
  return(parsed_obj)
}

#' Parse the Kimmerling et al. (2018) fibroblast data set
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @export
parse_Kimmerling <- function(absolute = TRUE, use_spike_ins = FALSE) {
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
  if(absolute) {
    for(i in 1:ncol(counts)) {
      # Multiply; library size should have a positive association with mass!
      counts[,i] <- counts[,i] * sf[i]
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
  
  counts <- as.matrix(counts)
  colnames(counts) <- NULL
  rownames(counts) <- NULL

  groups <- factor(groups, levels = c("1", "2"))
  levels(groups) <- c("low_mass", "high_mass")
  
  parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
  return(parsed_obj)
}

#' Parse the TCGA ESCA RNA-seq data
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @export
parse_ESCA <- function(absolute = TRUE) {
  file_dir <- file.path("data", "TCGA_ESCA")
  
  # Think about how to even utilize spike-ins for this normalization
  # See: https://www.biostars.org/p/81803/
  # Or Lun et al. 2017
  
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

  # Find spike-ins and separate into spike-in and all other sequence data sets
  # spikein_idx <- unname(which(sapply(esca$gene_name, function(x) {
  #   str_detect(x, "^ERCC")
  # })))
  esca_spike <- esca[spikein_idx,2:ncol(esca)]
  esca <- esca[-spikein_idx,2:ncol(esca)]
  
  # Label samples with tumor subtype
  barcodes <- unname(sapply(colnames(esca), function(x) {
    str_replace_all(x, "\\.", "-")
  }))
  scc_flag <- rep("Other", ncol(esca))
  scc_flag[barcodes %in% sccs$Patient.Barcodes] <- "SCC"
  
  # Visualizing via PCA show these two tumor types are very different
  # temp <- esca[sample(1:nrow(esca), size = 200),]
  # 
  # coords <- cmdscale(dist(t(temp)))
  # ggplot(data.frame(x = coords[,1], y = coords[,2], scc = scc_flag),
  #        aes(x = x, y = y, fill = factor(scc))) +
  #   geom_point(size = 2, shape = 21)
  
  # Plot omitted here but there's an overall higher spike-in abundance in the
  # SCC samples
  
  if(absolute) {
    size_factor <- unname(colMeans(esca_spike))
    counts <- esca / size_factor
    counts <- 2**counts - 1 # assume a pseudocount of 1 was added here
    groups <- scc_flag
    tax = NULL
  } else {
    counts <- 2**esca - 1 # assume a pseudocount of 1 was added here
    groups <- scc_flag
    tax = NULL
  }
  
  # Clean up all the labels we won't use
  rownames(counts) <- NULL
  colnames(counts) <- NULL
  counts <- as.matrix(counts)

  counts <- cbind(counts[,groups == "Other"], counts[,groups == "SCC"])
  groups <- sort(groups)
  
  parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
  
  # saveRDS(parsed_obj, file = parsed_data_fn)
  return(parsed_obj)
}


#' Parse the Vieira-Silva et al. (2019) QMP data set
#'
#' @param absolute flag indicating whether or not to parse absolute abundance
#' data
#' @return named list of counts and group labels
#' @import stringr
#' @export
parse_VieiraSilva <- function(absolute = TRUE) {
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
    left_join(data.frame(sample_id = rownames(data), idx_in_data = 1:nrow(data)), by = "sample_id")
  data <- data[mapping$idx_in_data,]
  data <- t(data)

  # Columns (samples) in `data` now match identity of rows in `meta`
  # Difference in total abundance is predictably largest between Crohn's disease
  # cohort (CD) and healthy controls (mHC)
  # library(ggplot)
  # ggplot(data.frame(y = colSums(counts), x = meta$cohort), aes(x = x, y = y)) +
  #   geom_boxplot()
  
  counts <- round(data)
  if(!absolute) {
    # Shuffle observed abundances
    set.seed(1001)
    new_totals <- sample(round(colSums(counts)))
    for(i in 1:ncol(counts)) {
      counts[,i] <- rmultinom(1, size = new_totals[i], prob = counts[,i] / sum(counts[,i]))
    }
  }
  counts <- as.matrix(counts)
  colnames(counts) <- NULL
  rownames(counts) <- NULL
  
  groups <- factor(meta$cohort)
  
  parsed_obj <- list(counts = counts, groups = groups, tax = NULL)
  return(parsed_obj)
}











