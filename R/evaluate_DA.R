#' Evaluate differential abundance with ANCOM-BC
#'
#' @param data simulated data set
#' @param groups group (cohort) labels
#' @param use_TMM if TRUE, uses the trimmed mean of M-values normalization
#' @return unadjusted p-values and "effect sizes" for DA for all features
#' @import phyloseq
#' @import ANCOMBC
#' @import stringr
#' @export
call_DA_ANCOMBC <- function(data, groups) {
  physeq <- phyloseq(otu_table(data, taxa_are_rows = FALSE),
                     sample_data(data.frame(group = groups)))
  out <- ancombc(phyloseq = physeq, formula = "group", p_adj_method = "BH")

  # Use log regression betas as an effect size in ANCOM-BC
  beta_df <- data.frame(feature = colnames(otu_table(physeq))) %>%
    left_join(data.frame(feature = rownames(out$res$beta),
                         beta = exp(unname(out$res$beta))), by = "feature")
  beta_df$beta[is.na(beta_df$beta)] <- 0
  
  # ANCOM will omit some very low abundance features from consideration
  # Re-incorporate these with p-values of 1
  pval_df <- data.frame(feature = colnames(otu_table(physeq))) %>%
    left_join(data.frame(feature = rownames(out$res$p_val),
                         pval = unname(out$res$p_val)), by = "feature")
  pval_df$pval[is.na(pval_df$pval)] <- 1
  
  # Combine these
  combo_df <- pval_df %>%
    left_join(beta_df, by = "feature")
  combo_df$feature <- str_replace(combo_df$feature, "sp", "gene")
  
  combo_df
}

#' Evaluate differential abundance with edgeR
#'
#' @param data simulated data set
#' @param groups group (cohort) labels
#' @param use_TMM if TRUE, uses the trimmed mean of M-values normalization
#' @return unadjusted p-values and "effect sizes" for DA for all features
#' @import edgeR
#' @export
call_DA_edgeR <- function(data, groups, use_TMM = FALSE) {
  # DGEList expects a features x samples orientation; accommodate this
  data <- t(data)
  n_genes <- nrow(data)
  # Estimate DE from relative abundances with edgeR and TMM
  dge_obj <- DGEList(counts = data, group = groups, lib.size = colSums(data))
  if(use_TMM) {
    dge_obj <- calcNormFactors(dge_obj) # optional bit
  }
  design <- model.matrix(~ groups)
  dge_obj <- estimateGLMCommonDisp(dge_obj, design)
  dge_obj <- estimateGLMTagwiseDisp(dge_obj, design)
  fit <- glmFit(dge_obj, design)
  lrt <- glmLRT(fit, coef = 2)
  pval_df <- data.frame(feature = paste0("gene", 1:nrow(data)),
             pval = lrt$table$PValue,
             beta = exp(unname(unlist(fit$coefficients[,2]))))
  pval_df
}

#' Evaluate differential abundance with a negative binomial GLM (via MASS)
#'
#' @param data simulated data set (samples x features)
#' @param groups group (cohort) labels
#' @return unadjusted p-values and "effect sizes" for DA for all features
#' @import MASS
#' @import dplyr
#' @export
call_DA_NB <- function(data, groups) {
  pval_df <- data.frame(feature = paste0("gene", 1:ncol(data)),
                        pval = 1,
                        beta = 0)
  data <- spike_in_ones(data, groups)
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
      pval_df[feature_idx,]$beta <- exp(coef(summary(fit))[2,1])
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
#' @param control_indices if specified, these features are used as the stable 
#' features against which to normalize observed abundances (DESeq2-only)
#' @return unadjusted p-values and "effect sizes" for DA for all features
#' @import DESeq2
#' @import dplyr
#' @export
call_DA_DESeq2 <- function(data, groups,
                           control_indices = NULL) {
  sampleData <- data.frame(group = factor(groups))
  dds <- suppressMessages(DESeqDataSetFromMatrix(countData = t(data), colData = sampleData, design = ~ group))
  if(!is.null(control_indices)) {
    dds <- suppressMessages(DESeq2::estimateSizeFactors(object = dds, controlGenes = control_indices))
  } else {
    dds <- suppressMessages(DESeq2::estimateSizeFactors(object = dds))
  }
  dds <- suppressMessages(DESeq2::estimateDispersions(object = dds, fitType = "local"))
  dds <- suppressMessages(DESeq2::nbinomWaldTest(object = dds))
  res <- DESeq2::results(object = dds, alpha = 0.05)
  pval_df <- data.frame(feature = paste0("gene", 1:ncol(data)),
                        pval = res$pvalue,
                        beta = 2**(res$log2FoldChange))
  pval_df$pval[is.na(pval_df$pval)] <- 1
  pval_df$beta[is.na(pval_df$beta)] <- 0
  pval_df
  
  # Previously
  # I had been using a Seurat wrapper when many more differential abundance
  # testing methods - including MAST - were in use. FindMarkers() was a 
  # conventient general purpose function for calling these.
  # count_table <- t(data)
  # n_genes <- nrow(count_table)
  # n_samples <- ncol(count_table)
  # 
  # # Create count table and metadata objects  
  # cell_metadata <- data.frame(active.ident = as.factor(groups))
  # rownames(cell_metadata) <- colnames(count_table)
  # rownames(count_table) <- paste0("gene", 1:n_genes)
  # colnames(count_table) <- paste0("cell", 1:n_samples)
  # levels(cell_metadata$active.ident) <- c("untreated", "treated")
  # rownames(cell_metadata) <- colnames(count_table)
  # 
  # # Create Seurat object manually:
  # # https://learn.gencore.bio.nyu.edu/single-cell-rnaseq/loading-your-own-data-in-seurat-reanalyze-a-different-dataset/
  # seurat_data <- CreateSeuratObject(counts = count_table,
  #                                   project = "simulation",
  #                                   min.cells = 1,
  #                                   min.features = 10,
  #                                   meta.data = cell_metadata)
  # DefaultAssay(seurat_data) <- "RNA"
  # # Assign labels manually
  # seurat_data <- SetIdent(seurat_data, value = seurat_data@meta.data$active.ident)
  # results <- suppressMessages(FindMarkers(seurat_data,
  #                                         ident.1 = "untreated",
  #                                         ident.2 = "treated",
  #                                         test.use = "DESeq2"))
  # pval_df <- left_join(data.frame(feature = rownames(count_table)),
  #                      data.frame(feature = rownames(results),
  #                                 pval = results$p_val),
  #                      by = "feature")
  # # It's likely some features will have been excluded from consideration if they
  # # are near-zero. Plug p-values of 1 in for these.
  # pval_df$pval[is.na(pval_df$pval)] <- 1
  # pval_df
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
#' @return unadjusted p-values and "effect sizes" for DA for all features
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
                        pval = results$we.ep,
                        beta = exp(results$effect))
  
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
#' @return unadjusted p-values and "effect sizes" for DA for all features
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
  
  # getMarkerEffects(results) doesn't seem to work - returns empty list
  # I'll do something a little cruder: scale the counts by hand (by scran's
  # computed scaling factor) and take the mean fold change across groups
  sf <- sizeFactors(sce)
  norm_counts <- counts(sce)
  for(i in 1:ncol(norm_counts)) {
    norm_counts[,i] <- norm_counts[,i]/sf[i]
  }
  groups <- factor(groups)
  est <- numeric(nrow(norm_counts))
  for(i in 1:nrow(norm_counts)) {
    est[i] <- mean(norm_counts[i,groups == levels(groups)[2]])/mean(norm_counts[i,groups == levels(groups)[1]])
  }
  pval_df$beta <- est
  
  pval_df
}

#' Wrapper for call_DA_ANCOMBC generating reference and target calls
#'
#' @param ref_data absolute abundance data set (samples x features)
#' @param data relative abundance data set (samples x features)
#' @param groups group (cohort) labels
#' @param oracle_calls optional baseline (true) differential abundance calls
#' @return baseline (true) differential abundance calls and calls made on 
#' relative abundances
#' @export
DA_by_ANCOMBC <- function(ref_data, data, groups, oracle_calls = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_NB(ref_data, groups)
  }
  calls <- call_DA_ANCOMBC(data, groups)
  return(list(oracle_calls = oracle_calls, calls = calls))
}

#' Wrapper for call_DA_edgeR generating reference and target calls
#'
#' @param ref_data absolute abundance data set (samples x features)
#' @param data relative abundance data set (samples x features)
#' @param groups group (cohort) labels
#' @param oracle_calls optional baseline (true) differential abundance calls
#' @param use_TMM if TRUE, uses the trimmed mean of M-values normalization
#' @return baseline (true) differential abundance calls and calls made on 
#' relative abundances
#' @export
DA_by_edgeR <- function(ref_data, data, groups, oracle_calls = NULL,
                        use_TMM = FALSE) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_NB(ref_data, groups)
  }
  calls <- call_DA_edgeR(data, groups, use_TMM = use_TMM)
  return(list(oracle_calls = oracle_calls, calls = calls))
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
#' @param control_indices if specified, these features are used as the stable 
#' features against which to normalize observed abundances (DESeq2-only)
#' @return baseline (true) differential abundance calls and calls made on 
#' relative abundances
#' @export
DA_by_DESeq2 <- function(ref_data, data, groups, oracle_calls = NULL,
                         control_indices = NULL) {
  if(is.null(oracle_calls)) {
    oracle_calls <- call_DA_DESeq2(ref_data, groups)
  }
  if(!is.null(control_indices)) {
    calls <- call_DA_DESeq2(data, groups,
                            control_indices = control_indices)
  } else {
    calls <- call_DA_DESeq2(data, groups)
  }
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
#' @param control_indices if specified, these features are used as the stable 
#' features against which to normalize observed abundances (DESeq2-only)
#' @return named list of baseline (true) differential abundance calls and
#' calls made on relative abundances (the observed data)
#' @export
DA_wrapper <- function(ref_data, data, groups, method = "NBGLM",
                       oracle_calls = NULL,
                       control_indices = NULL) {
  if(method == "NBGLM") {
    DA_calls <- DA_by_NB(ref_data, data, groups, oracle_calls = oracle_calls)
  }
  if(method == "ANCOMBC") {
    DA_calls <- DA_by_ANCOMBC(ref_data, data, groups, oracle_calls = oracle_calls)
  }
  if(method == "DESeq2") {
    DA_calls <- DA_by_DESeq2(ref_data, data, groups, oracle_calls = oracle_calls,
                             control_indices = control_indices)
  }
  if(method == "MAST") {
    DA_calls <- DA_by_MAST(ref_data, data, groups, oracle_calls = oracle_calls)
  }
  if(method == "ALDEx2") {
    DA_calls <- DA_by_ALDEx2(ref_data, data, groups,
                             oracle_calls = oracle_calls)
  }
  if(method == "scran") {
    DA_calls <- DA_by_scran(ref_data, data, groups, oracle_calls = oracle_calls)
  }
  if(method == "edgeR") {
    DA_calls <- DA_by_edgeR(ref_data, data, groups,
                            oracle_calls = oracle_calls, use_TMM = FALSE)
  }
  if(method == "edgeR_TMM") {
    DA_calls <- DA_by_edgeR(ref_data, data, groups,
                            oracle_calls = oracle_calls, use_TMM = TRUE)
  }
  if(is.null(oracle_calls)) {
    # oracle_calls <- DA_calls$oracle_calls$pval
    oracle_calls <- DA_calls$oracle_calls
  }
  # calls <- DA_calls$calls$pval
  calls <- DA_calls$calls
  
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
#' @param alpha significance threshold after FD adjustment
#' @param beta optional effect size threshold at the scale of the original 
#' abundances (i.e. not log scale)
#' @return named list of all calls, TPR, and FPR
#' @export
calc_DA_discrepancy <- function(calls, oracle_calls, adjusted = TRUE,
                                alpha = 0.05, beta = NULL) {
  if(typeof(calls$pval) == "character") {
    calls$pval <- as.numeric(strsplit(calls$pval, ";")[[1]])
  }
  if(typeof(calls$beta) == "character") {
    calls$beta <- as.numeric(strsplit(calls$beta, ";")[[1]])
  }
  if(typeof(oracle_calls$pval) == "character") {
    oracle_calls$pval <- as.numeric(strsplit(oracle_calls$pval, ";")[[1]])
  }
  if(typeof(oracle_calls$beta) == "character") {
    oracle_calls$beta <- as.numeric(strsplit(oracle_calls$beta, ";")[[1]])
  }
  
  if(adjusted) {
    calls$pval <- p.adjust(calls$pval, method = "BH")
    oracle_calls$pval <- p.adjust(oracle_calls$pval, method = "BH")
  }
  
  if(!is.null(beta)) {
    de <- calls$pval < alpha & (calls$beta < 1/beta | calls$beta > beta)
    sim_de <- oracle_calls$pval < alpha & (oracle_calls$beta < 1/beta | oracle_calls$beta > beta)
  } else {
    de <- calls$pval < alpha
    sim_de <- oracle_calls$pval < alpha
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
