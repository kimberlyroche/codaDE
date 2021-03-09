#' Apply scran normalization and call differentially expressed (marker) genes
#' 
#' @param data simulated data set
#' @param groups group (cohort) labels
#' @import SingleCellExperiment
#' @import scran
#' @import scater
#' @export
call_DA_scran <- function(data, groups) {
  # data <- sim_data$observed_counts1[,1:100]
  count_table <- t(data)
  # groups <- sim_data$groups
  n_genes <- nrow(count_table)
  n_samples_condition <- ncol(count_table)/2
  cell_metadata <- data.frame(condition = groups)
  
  rownames(cell_metadata) <- colnames(count_table)
  rownames(count_table) <- paste0("gene", 1:n_genes)
  colnames(count_table) <- paste0("cell", 1:(n_samples_condition*2))
  
  sce <- SingleCellExperiment(assays = list(counts = count_table),
                              colData = cell_metadata)

  # Note: We'll need to add QC stuff here ultimately.
  
  # Methods below taken from this tutorial:
  # https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html
  
  clusters <- suppressWarnings(quickCluster(sce, min.size = 1))

  sce <- suppressWarnings(computeSumFactors(sce, clusters = clusters)) # Use a first-pass clustering
  sce <- logNormCounts(sce)
  
  # Cluster from docs: https://rdrr.io/bioc/scran/man/getClusteredPCs.html
  sce <- suppressWarnings(scater::runPCA(sce))
  output <- getClusteredPCs(reducedDim(sce))
  
  npcs <- metadata(output)$chosen # number of "informative dimensions"
  reducedDim(sce, "PCAsub") <- reducedDim(sce, "PCA")[,1:npcs,drop=FALSE]
  
  # Build a graph based on this truncated dimensionality reduction
  g <- buildSNNGraph(sce, use.dimred = "PCAsub")
  cluster <- igraph::cluster_walktrap(g)$membership
  
  # Assigning to the 'colLabels' of the 'sce'
  colLabels(sce) <- factor(cluster)
  # table(colLabels(sce))
  
  # Visualize the clustering
  # sce <- runTSNE(sce, dimred = "PCAsub")
  # plotTSNE(sce, colour_by = "label", text_by = "label")
  
  markers <- findMarkers(sce, test.type = "wilcox") # t-test by default
  filter_vector <- markers[[1]][,1:3]$FDR < 0.05
  sig_genes <- rownames(markers[[1]])[filter_vector]
  calls <- logical(n_genes)
  calls[which(rownames(count_table) %in% sig_genes)] <- TRUE
  return(calls)
}

#' Fit negative binomial model to null and full models and evaluate differential abundance for a focal gene
#'
#' @param data simulated data set
#' @param groups group (cohort) labels
#' @param feature_idx index of gene to test for differential abundance
#' @return p-value from NB GLM fit with MASS::glm.nb associated with group coefficient
#' @import MASS
#' @export
call_DA_NB <- function(data, groups, feature_idx) {
  gene_data <- data.frame(counts = data[,feature_idx], groups = groups)
  fit <- tryCatch({
    glm.nb(counts ~ groups, data = gene_data)
  }, warning = function(w) {
    # this is typically "iteration limit reached"
  }, error = function(err) {
    # this is typically "NaNs produced", produced from under/overflow and ultimately Poisson-like dispersion
    # https://r.789695.n4.nabble.com/R-Error-Warning-Messages-with-library-MASS-using-glm-td4556771.html
    # using a quasipoisson here seems to work well
  })
  if(is.null(fit)) {
    fit <- tryCatch({
      glm(counts ~ groups, family = quasipoisson(), data = gene_data)
    }, warning = function(w) {}, error = function(err) { })
  }
  if(!is.null(fit)) {
    return(coef(summary(fit))[2,4])
  } else {
    cat(paste0("Fit failed on feature",feature_idx,"\n"))
    return(NA)
  }
}

#' Fit log linear model to null and full models and evaluate differential abundance for a focal gene
#'
#' @param data simulated data set
#' @param feature_idx index of gene to test for differential abundance
#' @param call_abundances if TRUE, call DA on abundances; if FALSE, on observed counts
#' @return p-value associated with group coefficient
#' @import lmPerm
#' @export
call_DA_LM <- function(data, feature_idx, call_abundances = TRUE) {
  if(call_abundances) {
    gene_data <- data.frame(log_counts = log(data$abundances[,feature_idx] + 0.5), groups = data$groups)
  } else {
    gene_data <- data.frame(log_counts = log(data$observed_counts[,feature_idx] + 0.5), groups = data$groups)
  }
  # there seems to be no way to suppress the status output
  discard_str <- capture.output(fit <- lmp(log_counts ~ groups, data = gene_data))
  pval <- summary(fit)$coef[2,3]
  return(pval)
}

#' Evaluate differential abundance with edgeR
#' This evaluates expression on all features of the count matrix together
#'
#' @param data simulated data set
#' @param call_abundances if TRUE, call DA on abundances; if FALSE, on observed counts
#' @param normalization_method if NULL, uses library size normalization; other options include "TMM" and "scran"
#' @return p-value associated with group coefficient
#' @import edgeR
#' @import SingleCellExperiment
#' @import scran
#' @export
call_DA_edgeR <- function(data, call_abundances = TRUE, normalization_method = NULL) {
  # DGEList expects samples as columns
  if(!is.null(normalization_method)) {
    if(normalization_method == "TMM") {
      if(call_abundances) {
        dge_obj <- DGEList(counts = t(data$abundances), group = factor(data$groups))
      } else {
        dge_obj <- DGEList(counts = t(data$observed_counts), group = factor(data$groups))
      }
      dge_obj <- calcNormFactors(dge_obj, method = "TMM")
    } else if(normalization_method == "scran") {
      # First, get the data into a SingleCellExperiment (Bioconductor) object
      if(call_abundances) {
        y <- t(data$abundances)
      } else {
        y <- t(data$observed_counts)
      }
      # These labels may not be strictly necessary but it's probably not bad form
      rownames(y) <- paste0("gene_", 1:nrow(y))
      y <- as.matrix(y)
      cell_metadata <- data.frame(condition = data$groups)
      rownames(cell_metadata) <- paste0("cell_", 1:ncol(y))
      sce <- SingleCellExperiment(assays = list(counts = y),
                                  colData = cell_metadata)
      clust.sce <- quickCluster(sce)
      sce <- suppressWarnings(scran::computeSumFactors(sce, cluster = clust.sce))
      dge_obj <- convertTo(sce, type = "edgeR")
    } else {
      stop("Unknown normalization method!")
    }
  } else {
    if(call_abundances) {
      dge_obj <- DGEList(counts = t(data$abundances), group = factor(data$groups))
    } else {
      dge_obj <- DGEList(counts = t(data$observed_counts), group = factor(data$groups))
    }
    dge_obj <- calcNormFactors(dge_obj, method = "none") # library size normalization-only
  }
  design <- model.matrix(~ data$groups)
  dge_obj <- estimateDisp(dge_obj, design)
  # LRT recommended for single-cell data
  fit <- glmFit(dge_obj, design)
  lrt <- glmLRT(fit, coef = 2)
  # alternatively: is.de <- decideTestsDGE(lrt)
  pval <- lrt@.Data[[14]]$PValue
  # quasi-likelihood recommended for bulk RNA-seq (different dispersion estimation procedure)
  # fit <- glmQLFit(dge_obj, design)
  # lrt <- glmQLFTest(fit, coef=2)
  # pval <- lrt@.Data[[17]]$PValue
  return(pval)
}

#' Get additive logratios from the data set
#'
#' @param data simulated data set
#' @param call_abundances if TRUE, call DA on abundances; if FALSE, on observed counts
#' @return logratios
#' @import driver
#' @export
get_logratios <- function(data, call_abundances = TRUE) {
  if(call_abundances) {
    # does it make sense to call these with NB or ALR-normal model?
    counts <- data$abundances
  } else {
    # use median abundance feature as the ALR reference
    counts <- data$observed_counts
  }
  logratios <- driver::alr(counts + 0.5)
  return(logratios)
}

#' Fit normal model to logratios and evaluate differential abundance for a focal gene
#' NOTE: probably need something better than a normal here; it's not a great fit :(
#'       predict from fitted model for a few of these and see if this is the case
#'
#' @param logratios data in ALR format
#' @param feature_idx index of gene to test for differential abundance
#' @param groups treatment group assignments across samples
#' @return p-value from lm
#' @export
call_DA_CODA <- function(logratios, feature_idx, groups) {
  # adjust the indexing
  gene_data <- data.frame(logratios = logratios[,feature_idx], groups = groups)
  fit <- lm(logratios ~ groups, data = gene_data)
  pval <- coef(summary(fit))[2,4]
  return(pval)
}

#' Generate simulated data with RNA-seq like features and > 3 fold differential abundance, evaluate error in
#' differential abundance calls and write these to an output file
#'
#' @param data simulated data set
#' @param alpha significant level below which to call a feature differentially abundant
#' @param use_ALR if TRUE, evaluates differential abundance using a spike-in and the additive logratio
#' @param filter_abundance the minimum average abundance to require of features we'll evaluate for differential abundance
#' (e.g. a filter_abundance of 1 evaluates features with average abundance across conditions of at least 1)
#' @param method differential abundance testing method to use; options are "NB", "GLM", "edgeR", "edgeR_TMM", "edgeR_scran"
#' @details Writes out error and simulation run statistics to a file.
#' @return NULL
#' @import edgeR
#' @export
evaluate_DA <- function(data, alpha = 0.05, use_ALR = FALSE, filter_abundance = 0, method = "edgeR") {
  # calculate and compare differential abundance calls
  cat("Evaluating differential abundance...\n")
  evaluate_features <- apply(data$observed_counts, 2, function(x) mean(x) > filter_abundance)
  p <- ncol(data$abundances)
  # Use intended DE as a baseline...
  calls.abundances <- rep(FALSE, p)
  calls.abundances[data$da_genes] <- TRUE
  calls.observed_counts <- rep(FALSE, p)
  
  if(use_ALR) {
    # this setup is currently only possible with NB DA calling
    logratios.abundances <- get_logratios(data, call_abundances = TRUE)
    logratios.observed_counts <- get_logratios(data, call_abundances = FALSE)
    for(i in 1:p) {
      if(evaluate_features[i]) {
        # pval.abundances <- call_DA_CODA(logratios.abundances, i, data$groups)
        pval.observed_counts <- call_DA_CODA(logratios.observed_counts, i, data$groups)
        #if(!is.na(pval.abundances) & !is.na(pval.observed_counts)) {
        if(!is.na(pval.observed_counts)) {
          # calls.abundances <- c(calls.abundances, pval.abundances <= alpha/length(evaluate_features))
          calls.observed_counts <- c(calls.observed_counts, pval.observed_counts <= alpha/length(evaluate_features))
        }
      }
    }
  } else {
    if(method %in% c("edgeR", "edgeR_TMM", "edgeR_scran")) {
      # Note: This inherently looks for DE over relative abundances!
      # pval.abundances <- call_DA_edgeR(data, call_abundances = TRUE)
      # calls.abundances <- pval.abundances <= alpha/length(evaluate_features)
      if(method == "edgeR_TMM") {
        pval.observed_counts <- call_DA_edgeR(data, call_abundances = FALSE, normalization_method = "TMM")
      } else if(method == "edgeR_scran") {
        pval.observed_counts <- call_DA_edgeR(data, call_abundances = FALSE, normalization_method = "scran")
      } else {
        pval.observed_counts <- call_DA_edgeR(data, call_abundances = FALSE)
      }
      calls.observed_counts <- pval.observed_counts <= alpha/length(evaluate_features)
    } else {
      for(i in 1:p) {
        if(i %% 100 == 0) {
          cat("Evaluating DA on feature:",i,"\n")
        }
        if(evaluate_features[i]) {
          if(method == "NB") {
            # pval.abundances <- call_DA_NB(data, i, call_abundances = TRUE)
            pval.observed_counts <- call_DA_NB(data, i, call_abundances = FALSE)
          } else if(method == "GLM") {
            # pval.abundances <- call_DA_LM(data, i, call_abundances = TRUE)
            pval.observed_counts <- call_DA_LM(data, i, call_abundances = FALSE)
          } else {
            stop("Unknown differential abundance calling method!")
          }
          #if(!is.na(pval.abundances) & !is.na(pval.observed_counts)) {
          if(!is.na(pval.observed_counts)) {
            # calls.abundances <- c(calls.abundances, pval.abundances <= alpha/length(evaluate_features))
            # calls.observed_counts <- c(calls.observed_counts, pval.observed_counts <= alpha/length(evaluate_features))
            calls.observed_counts[i] <- pval.observed_counts <= alpha/length(evaluate_features)
          }
        }
      }
    }
  }
  
  FN <- sum(calls.abundances == TRUE & calls.observed_counts == FALSE)
  FP <- sum(calls.abundances == FALSE & calls.observed_counts == TRUE)
  TN <- sum(calls.abundances == FALSE & calls.observed_counts == FALSE)
  TP <- sum(calls.abundances == TRUE & calls.observed_counts == TRUE)
  quantity_evaluated <- "counts"
  if(use_ALR) {
    quantity_evaluated <- "logratios"
  }
  
  library_size.abundances = rowSums(data$abundances)
  library_size.observed_counts = rowSums(data$observed_counts)
  
  return(list(FN = FN,
              FP = FP,
              TN = TN,
              TP = TP,
              TPR = TP / (TP + FN),
              FPR = FP / (FP + TN),
              quantity_evaluated = quantity_evaluated,
              method = method,
              filter_abundance = filter_abundance,
              no_features_above_threshold = sum(evaluate_features),
              no_features_detectable = sum(calls.abundances == TRUE),
              calls.abundances = calls.abundances,
              calls.observed_counts = calls.observed_counts,
              library_size.abundances = library_size.abundances,
              library_size.observed_counts = library_size.observed_counts))
}

#' Generate simulated data with RNA-seq like features and > 3 fold differential abundance, evaluate error in
#' differential abundance calls and write these to an output file
#'
#' @param p number of genes
#' @param n number of samples per group
#' @param proportion_da proportion of differentially abundant genes
#' @param k number of cell types to simulate
#' @param library_size_correlation correlation of observed abundances to original total abundances
#' @param alpha significant level below which to call a feature differentially abundant
#' @param use_ALR if TRUE, evaluates differential abundance using a spike-in and the additive logratio
#' @param filter_abundance the minimum average abundance to require of features we'll evaluate for differential abundance
#' (e.g. a filter_abundance of 1 evaluates features with average abundance across conditions of at least 1)
#' @param methods vector of differential abundance testing methods to use; options are "NB", "GLM", "edgeR"
#' @details Writes out error and simulation run statistics to a file.
#' @return NULL
#' @import entropy
#' @import uuid
#' @export
run_RNAseq_evaluation_instance <- function(p, n, proportion_da, k = 1, library_size_correlation = 0,
                                           alpha = 0.05, use_ALR = FALSE, filter_abundance = 1,
                                           methods = c("NB", "edgeR")) {
  data <- simulate_sequence_counts(p = p, n = n, k = k, proportion_da = proportion_da, library_size_correlation = library_size_correlation,
                                     spike_in = use_ALR)
  # we've got to make sure there's non-zero variation in all genes in each condition
  # usually this fails to be true if one gene is 100% unobserved in one condition and minimally present in the other
  # as a small workaround, if we find any genes will fully no expression in a given condition, we'll drop in a
  #   single 1-count
  # this shouldn't change results and will allow us to fit the GLM across all genes
  
  # which genes are fully missing in the original counts for group 0?
  drop_in <- which(colSums(data$abundances[data$groups == 0,]) == 0)
  for(d in drop_in) {
    data$abundances[sample(which(data$groups == 0))[1],d] <- 1
  }
  
  # which genes are fully missing in the original counts for group 0?
  drop_in <- which(colSums(data$abundances[data$groups == 1,]) == 0)
  for(d in drop_in) {
    data$abundances[sample(which(data$groups == 1))[1],d] <- 1
  }
  
  # which genes are fully missing in the original counts for group 0?
  drop_in <- which(colSums(data$observed_counts[data$groups == 0,]) == 0)
  for(d in drop_in) {
    data$observed_counts[sample(which(data$groups == 0))[1],d] <- 1
  }
  
  # which genes are fully missing in the original counts for group 0?
  drop_in <- which(colSums(data$observed_counts[data$groups == 0,]) == 0)
  for(d in drop_in) {
    data$observed_counts[sample(which(data$groups == 0))[1],d] <- 1
  }
  
  for(method in methods) {
    save_obj <- evaluate_DA(data, alpha = alpha, use_ALR = use_ALR, filter_abundance = filter_abundance, method = method)
    data_obj <- save_obj
    data_obj$p <- p
    data_obj$proportion_da <- proportion_da
    data_obj$library_size_correlation <- library_size_correlation
    new_filename <- paste0(UUIDgenerate(), ".rds")
    # saveRDS(data_obj, file.path("simulated_data", new_filename))
    save_path <- file.path("simulated_analyses", new_filename)
    cat("Saving to",save_path,"\n")
    saveRDS(data_obj, save_path)
  }
}

#' Generate simulated data with RNA-seq like features and > 3 fold differential abundance, evaluate false negatives and positives in
#' differential abundance calls and write these to an output file
#'
#' @param p number of genes
#' @param n number of samples per group
#' @param k number of cell types to simulate
#' @param de_sweep vector of proportions of differentially abundant genes to simulate
#' @param corr_sweep vector of correlations of size factors (true vs. observed abundances) to simulate
#' @param alpha significant level below which to call a feature differentially abundant
#' @param use_ALR if TRUE, evaluates differential abundance using a spike-in and the additive logratio
#' @param filter_abundance the minimum average abundance to require of features we'll evaluate for differential abundance
#' (e.g. a filter_abundance of 1 evaluates features with average abundance across conditions of at least 1)
#' @param methods vector of differential abundance testing methods to use; options are "NB", "GLM", "edgeR"
#' @details Writes out error and simulation run statistics to a file.
#' @return NULL
#' @export
sweep_simulations <- function(p, n, k = 1, de_sweep = seq(from = 0.1, to = 0.9, by = 0.1),
                              corr_sweep = seq(from = 0.1, to = 0.9, by = 0.1), alpha = 0.05, use_ALR = FALSE,
                              filter_abundance = 0, methods = c("NB", "edgeR")) {
  for(de_prop in de_sweep) {
    for(sf_corr in corr_sweep) {
      out_str <- paste0("w/ ",p," genes, DA proportion = ",round(de_prop, 2),", and size factor correlation = ",round(sf_corr, 2),"\n")
      cat("Single-cell RNA-seq",out_str)
      run_RNAseq_evaluation_instance(p, n, proportion_da = de_prop, k = k,
                                     library_size_correlation = sf_corr,
                                     alpha = 0.05, use_ALR = use_ALR, filter_abundance = filter_abundance,
                                     methods = methods)
    }
  }
}
