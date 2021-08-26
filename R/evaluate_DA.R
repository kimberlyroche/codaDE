#' #' Evaluate differential abundance with edgeR
#' #' This evaluates expression on all features of the count matrix together
#' #'
#' #' @param data simulated data set
#' #' @param groups group (cohort) labels
#' #' @param normalization_method if NULL, uses library size normalization; 
#' other options include "TMM" and "scran"
#' #' @return p-value for DA for all features
#' #' @import edgeR
#' #' @export
#' call_DA_edgeR <- function(data, groups, normalization_method = NULL) {
#'   # DGEList expects samples as columns
#'   if(!is.null(normalization_method)) {
#'     if(normalization_method == "TMM") {
#'       dge_obj <- DGEList(counts = t(data), group = factor(groups))
#'       dge_obj <- calcNormFactors(dge_obj, method = "TMM")
#'     } else {
#'       stop("Unknown normalization method!")
#'     }
#'   } else {
#'     dge_obj <- DGEList(counts = t(data), group = factor(groups))
#'     dge_obj <- calcNormFactors(dge_obj, method = "none") # library size 
#'     normalization-only
#'   }
#'   design <- model.matrix(~ groups)
#'   dge_obj <- estimateDisp(dge_obj, design)
#'   # LRT recommended for single-cell data
#'   fit <- glmFit(dge_obj, design)
#'   lrt <- glmLRT(fit, coef = 2)
#'   # alternatively: is.de <- decideTestsDGE(lrt)
#'   pval <- lrt@.Data[[14]]$PValue
#'   # quasi-likelihood recommended for bulk RNA-seq (different dispersion 
#'   estimation procedure)
#'   # fit <- glmQLFit(dge_obj, design)
#'   # lrt <- glmQLFTest(fit, coef=2)
#'   # pval <- lrt@.Data[[17]]$PValue
#'   return(pval)
#' }

#' Evaluate differential abundance with a negative binomial GLM (via MASS)
#'
#' @param data simulated data set (samples x features)
#' @param groups group (cohort) labels
#' @return unadjusted p-values for all features
#' @import MASS
#' @import dplyr
#' @export
call_DA_NB <- function(data, groups) {
  pval_df <- data.frame(feature = paste0("feature", 1:ncol(data)),
                        pval = 1)
  for(feature_idx in 1:ncol(data)) {
    gene_data <- data.frame(counts = data[,feature_idx], groups = groups)
    fit <- tryCatch({
      glm.nb(counts ~ groups, data = gene_data)
    }, warning = function(w) {
      # This is typically "iteration limit reached"
    }, error = function(err) {
      # This is typically "NaNs produced", produced from under/overflow and 
      # ultimately Poisson-like dispersion
      # https://r.789695.n4.nabble.com/R-Error-Warning-Messages-with-library-MASS-using-glm-td4556771.html
      # Using a quasipoisson here seems to work well
    })
    if(is.null(fit)) {
      fit <- tryCatch({
        glm(counts ~ groups, family = quasipoisson(), data = gene_data)
      }, warning = function(w) {}, error = function(err) { })
    }
    if(!is.null(fit)) {
      pval_df[feature_idx,]$pval <- coef(summary(fit))[2,4]
    } else {
      cat(paste0("Fit failed on feature",feature_idx,"\n"))
    }
  }
  pval_df
}

#' Evaluate differential abundance with DESeq2 (via Seurat)
#'
#' @param data simulated data set (samples x features)
#' @param groups group (cohort) labels
#' @return unadjusted p-values for all features
#' @import Seurat
#' @import dplyr
#' @export
call_DA_DESeq2 <- function(data, groups) {
  count_table <- t(data)
  n_genes <- nrow(count_table)
  n_samples <- ncol(count_table)
  
  # Create count table and metadata objects  
  cell_metadata <- data.frame(active.ident = as.factor(groups))
  rownames(cell_metadata) <- colnames(count_table)
  rownames(count_table) <- paste0("gene", 1:n_genes)
  colnames(count_table) <- paste0("cell", 1:n_samples)
  levels(cell_metadata$active.ident) <- c("untreated", "treated")
  rownames(cell_metadata) <- colnames(count_table)
  
  # Create Seurat object manually:
  # https://learn.gencore.bio.nyu.edu/single-cell-rnaseq/loading-your-own-data-in-seurat-reanalyze-a-different-dataset/
  seurat_data <- CreateSeuratObject(counts = count_table,
                                    project = "simulation",
                                    min.cells = 1,
                                    min.features = 10,
                                    meta.data = cell_metadata)
  DefaultAssay(seurat_data) <- "RNA"
  # Assign labels manually
  seurat_data <- SetIdent(seurat_data, value = seurat_data@meta.data$active.ident)
  results <- suppressMessages(FindMarkers(seurat_data,
                                          ident.1 = "untreated",
                                          ident.2 = "treated",
                                          test.use = "DESeq2"))
  pval_df <- left_join(data.frame(feature = rownames(count_table)),
                       data.frame(feature = rownames(results),
                                  pval = results$p_val),
                       by = "feature")
  # It's likely some features will have been excluded from consideration if they
  # are near-zero. Plug p-values of 1 in for these.
  pval_df$pval[is.na(pval_df$pval)] <- 1
  pval_df
}

#' Evaluate differential abundance with MAST (via Seurat)
#'
#' @param data simulated data set (samples x features)
#' @param groups group (cohort) labels
#' @return unadjusted p-values for all features
#' @import Seurat
#' @import dplyr
#' @export
call_DA_MAST <- function(data, groups) {
  count_table <- t(data)
  n_genes <- nrow(count_table)
  n_samples <- ncol(count_table)

  # Create count table and metadata objects  
  cell_metadata <- data.frame(active.ident = as.factor(groups))
  rownames(cell_metadata) <- colnames(count_table)
  rownames(count_table) <- paste0("gene", 1:n_genes)
  colnames(count_table) <- paste0("cell", 1:n_samples)
  levels(cell_metadata$active.ident) <- c("untreated", "treated")
  rownames(cell_metadata) <- colnames(count_table)
  
  # Create Seurat object manually, as in call_DA_DESeq2()
  seurat_data <- CreateSeuratObject(counts = count_table,
                                    project = "simulation",
                                    min.cells = 1,
                                    min.features = 10,
                                    meta.data = cell_metadata)
  DefaultAssay(seurat_data) <- "RNA"
  # Assign labels manually
  seurat_data <- SetIdent(seurat_data,
                          value = seurat_data@meta.data$active.ident)
  results <- suppressMessages(FindMarkers(seurat_data,
                                          ident.1 = "untreated",
                                          ident.2 = "treated",
                                          test.use = "MAST"))
  pval_df <- left_join(data.frame(feature = rownames(count_table)),
                       data.frame(feature = rownames(results),
                                  pval = results$p_val),
                       by = "feature")
  # It's likely some features will have been excluded from consideration if they
  # are near-zero. Plug p-values of 1 in for these.
  pval_df$pval[is.na(pval_df$pval)] <- 1
  pval_df
}

#' Evaluate differential abundance with ALDEx2
#'
#' @param data simulated data set (samples x features)
#' @param groups group (cohort) labels
#' @return unadjusted p-values for all features
#' @import ALDEx2
#' @import dplyr
#' @export
call_DA_ALDEx2 <- function(data, groups) {
  # Test using Welch's t-test (w/ IQLR representation)
  count_table <- t(data)
  n_genes <- nrow(count_table)
  n_samples <- ncol(count_table)
  rownames(count_table) <- paste0("gene", 1:n_genes)
  colnames(count_table) <- paste0("cell", 1:n_samples)
  
  results <- suppressMessages(aldex(count_table, groups, denom = "iqlr"))
  pval_df <- data.frame(feature = rownames(results),
                        pval = results$we.ep)
  
  # Interpreting the output:
  # we.ep - Expected P value of Welch's t test
  # we.eBH - Expected Benjamini-Hochberg corrected P value of Welch's t test
  # wi.ep - Expected P value of Wilcoxon rank test
  # wi.eBH - Expected Benjamini-Hochberg corrected P value of Wilcoxon test
  # kw.ep - Expected P value of Kruskal-Wallace test
  # kw.eBH - Expected Benjamini-Hochberg corrected P value of KW test
  
  pval_df
}

#' Evaluate differential abundance with scran
#'
#' @param data simulated data set (samples x features)
#' @param groups group (cohort) labels
#' @return unadjusted p-values for all features
#' @import SingleCellExperiment
#' @import scran
#' @import scater
#' @importFrom S4Vectors metadata
#' @export
call_DA_scran <- function(data, groups) {
  count_table <- t(data)
  n_genes <- nrow(count_table)
  cell_metadata <- data.frame(condition = groups)
  rownames(cell_metadata) <- colnames(count_table)
  rownames(count_table) <- paste0("gene", 1:n_genes)
  colnames(count_table) <- paste0("cell", 1:ncol(count_table))
  
  sce <- SingleCellExperiment(assays = list(counts = count_table),
                              colData = cell_metadata)
  
  # Methods below taken from this tutorial:
  # https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html
  # Use a first-pass clustering

  if(FALSE) {
    # (1) Find clusters
    clusters <- suppressWarnings(quickCluster(sce, min.size = 1))
  } else {
    # (2) Or specify them
    clusters <- as.numeric(as.factor(groups))
    if(min(clusters) == 0) {
      clusters <- clusters + 1
    }
  }
  sce <- suppressWarnings(computeSumFactors(sce, clusters = clusters))
  sce <- logNormCounts(sce)
  
  # ----------------------------------------------------------------------------
  #   Find clusters (1) or specify them (2)
  # ----------------------------------------------------------------------------
  if(FALSE) {
    # Cluster from docs: https://rdrr.io/bioc/scran/man/getClusteredPCs.html
    sce <- suppressWarnings(scater::runPCA(sce))
    output <- suppressWarnings(getClusteredPCs(reducedDim(sce)))
    npcs <- metadata(output)$chosen # number of "informative dimensions"
    reducedDim(sce, "PCAsub") <- reducedDim(sce, "PCA")[, 1:npcs, drop = FALSE]
    
    # Build a graph based on this truncated dimensionality reduction
    g <- suppressWarnings(buildSNNGraph(sce, use.dimred = "PCAsub"))
    cluster <- igraph::cluster_walktrap(g)$membership
    
    # Visualize the clustering
    # sce <- runTSNE(sce, dimred = "PCAsub")
    # plotTSNE(sce, colour_by = "label", text_by = "label")
  } else {
    clusters <- as.numeric(as.factor(groups))
    if(min(clusters) == 0) {
      clusters <- clusters + 1
    }
  }
  # Assigning to the 'colLabels' of the 'sce'
  colLabels(sce) <- factor(clusters)
  # table(colLabels(sce))
  
  # The FDR column contains BH-corrected p-values, which you can confirm by
  # calling plot(p.adjust(p.value, method = "BH", n = length(p.value)).
  markers <- findMarkers(sce, test.type = "wilcox") # t-test by default
  results <- markers[[1]]
  
  pval_df <- left_join(data.frame(feature = rownames(count_table)),
                       data.frame(feature = rownames(results),
                                  pval = results$p.value),
                       by = "feature")
  pval_df
}

#' Wrapper for call_DA_NB generating reference and target calls
#'
#' @param ref_data absolute abundance data set (samples x features)
#' @param data relative abundance data set (samples x features)
#' @param groups group (cohort) labels
#' @param oracle_calls optional baseline (true) differential abundance calls
#' @return baseline (true) differential abundance calls and calls made on 
#' relative abundances
#' @export
DA_by_NB <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_NB(ref_data, groups)
  }
  calls <- call_DA_NB(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

#' Wrapper for call_DA_DESeq2 generating reference and target calls
#'
#' @param ref_data absolute abundance data set (samples x features)
#' @param data relative abundance data set (samples x features)
#' @param groups group (cohort) labels
#' @param oracle_calls optional baseline (true) differential abundance calls
#' @return baseline (true) differential abundance calls and calls made on 
#' relative abundances
#' @export
DA_by_DESeq2 <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_DESeq2(ref_data, groups)
  }
  calls <- call_DA_DESeq2(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

#' Wrapper for call_DA_MAST generating reference and target calls
#'
#' @param ref_data absolute abundance data set (samples x features)
#' @param data relative abundance data set (samples x features)
#' @param groups group (cohort) labels
#' @param oracle_calls optional baseline (true) differential abundance calls
#' @return baseline (true) differential abundance calls and calls made on 
#' relative abundances
#' @export
DA_by_MAST <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_MAST(ref_data, groups)
  }
  calls <- call_DA_MAST(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

#' Wrapper for call_DA_ALDEx2 generating reference and target calls
#'
#' @param ref_data absolute abundance data set (samples x features)
#' @param data relative abundance data set (samples x features)
#' @param groups group (cohort) labels
#' @param oracle_calls optional baseline (true) differential abundance calls
#' @return baseline (true) differential abundance calls and calls made on 
#' relative abundances
#' @export
DA_by_ALDEx2 <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_ALDEx2(ref_data, groups)
  }
  calls <- call_DA_ALDEx2(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

#' Wrapper for call_DA_scran generating reference and target calls
#'
#' @param ref_data absolute abundance data set (samples x features)
#' @param data relative abundance data set (samples x features)
#' @param groups group (cohort) labels
#' @param oracle_calls optional baseline (true) differential abundance calls
#' @return baseline (true) differential abundance calls and calls made on 
#' relative abundances
#' @export
DA_by_scran <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_scran(ref_data, groups)
  }
  calls <- call_DA_scran(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

#' Wrapper function for differential abundance calling
#'
#' @param ref_data absolute abundance data set (samples x features)
#' @param data relative abundance data set (samples x features)
#' @param groups group (cohort) labels
#' @param method differential abundance calling method (e.g. "DESeq2")
#' @param oracle_calls optional baseline (true) differential abundance calls
#' @return named list of baseline (true) differential abundance calls and
#' calls made on relative abundances (the observed data)
#' @export
DA_wrapper <- function(ref_data, data, groups, method = "NBGLM",
                                oracle_calls = NULL) {
  if(method == "NBGLM") {
    DE_calls <- DA_by_NB(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
  if(method == "DESeq2") {
    DE_calls <- DA_by_DESeq2(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
  if(method == "MAST") {
    DE_calls <- DA_by_MAST(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
  if(method == "ALDEx2") {
    DE_calls <- DA_by_ALDEx2(ref_data, data, groups,
                             oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
  if(method == "scran") {
    DE_calls <- DA_by_scran(ref_data, data, groups, oracle_calls = oracle_calls)
    if(is.null(oracle_calls)) {
      oracle_calls <- DE_calls$oracle_calls$pval
    }
    calls <- DE_calls$calls$pval
  }
    
  return(list(calls = calls, oracle_calls = oracle_calls))
}

#' Calculates the discrepancy between differential abundance calls made on
#' absolute and relative abundances
#'
#' @param calls calls (as a numeric vector of p-values or a character string 
#' that can be parsed into one) made on relative abundances (the observed data)
#' @param oracle_calls calls (as a numeric vector of p-values or a character 
#' string that can be parsed into one) made on the absolute abundances
#' @param adjusted flag indicating whether or not to perform multiple test
#' correction
#' @return named list of all calls, TPR, and FPR
#' @export
calc_DA_discrepancy <- function(calls, oracle_calls, adjusted = TRUE) {
  if(typeof(calls) == "character") {
    calls <- as.numeric(strsplit(calls, ";")[[1]])
  }
  if(typeof(oracle_calls) == "character") {
    oracle_calls <- as.numeric(strsplit(oracle_calls, ";")[[1]])
  }
  
  if(adjusted) {
    # Calls on observed
    de <- p.adjust(calls, method = "BH") < 0.05
    # Calls on original abundances
    sim_de <- p.adjust(oracle_calls, method = "BH") < 0.05
  } else {
    de <- calls < 0.05
    sim_de <- oracle_calls < 0.05
  }

  TP_calls <- de & sim_de
  TP <- sum(TP_calls)
  FP_calls <- de & !sim_de
  FP <- sum(FP_calls)
  TN_calls <- !de & !sim_de
  TN <- sum(TN_calls)
  FN_calls <- !de & sim_de
  FN <- sum(FN_calls)
  
  # Occasionally all hits are spurious! Prevent NaN's here by interpreting the
  # TPR as zero when this happens.
  if(TP == 0 && FN == 0) {
    TPR <- 0
  } else {
    TPR <- TP/(TP+FN)
  }
  FPR <- FP/(FP+TN)

  # Special case: if there is no detectable differential abundance in the
  # reference, exclude this case altogether.
  if(sum(sim_de) == 0) {
    TPR <- NA
    FPR <- NA
  }
  
  return(list(TP_calls = TP_calls,
              FP_calls = FP_calls,
              TN_calls = TN_calls,
              FN_calls = FN_calls,
              TPR = TPR,
              FPR = FPR))
}

#' Threshold to generate calls on differentially abundant features
#'
#' @param counts abundance data set (samples x features)
#' @param fc_lower a lower threshold on fold change; smaller than this and 
#' observed change will be recorded as "differential"
#' @param fc_upper an upper threshold on fold change; larger than this and
#' observed change will be recorded as "differential"
#' @param nA optional parameter specifying the number of samples in condition A
#' @return differential calls in the form of faux p-values (0 for differential,
#' 1 for not differential)
#' @export
calc_threshold_DA <- function(counts, fc_lower = 0.5, fc_upper = 1.5, nA = NULL) {
  if(is.null(nA)) {
    nA <- nrow(counts)/2
    nB <- nA
  }
  m1 <- colMeans(counts[1:nA,] + 0.1)
  m2 <- colMeans(counts[(nA+1):nrow(counts),] + 0.1)
  fc <- m1 / m2
  oracle_calls <- rep(1, ncol(counts))
  oracle_calls[fc <= fc_lower | fc >= fc_upper] <- 0
  oracle_calls
}















