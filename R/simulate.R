# utility function to simulate a vector with some desired correlation to a reference vector
# source: https://stackoverflow.com/questions/47775389/how-to-simulate-a-vector-that-is-correlated-in-a-different-way-to-two-other-ex
simcor <- function(x, ymean = 0, ysd = 1, correlation = 0) {
  n <- length(x)
  y <- rnorm(n)
  z <- correlation * scale(x)[,1] + sqrt(1 - correlation^2) * scale(resid(lm(y ~ x)))[,1]
  yresult <- ymean + ysd * z
  yresult
}

#' Simulate sequence counts according to a negative binomial
#'
#' @param theta vector of (2) group dispersion parameters and (2) per-gene log covariates is theta
#' @param groups condition (e.g. control vs. treatment) labels
#' @param size_factor_correlation correlation of observed abundances to original absolute abundances
#' @param spike_in flag indicating whether or not to simulate a spike-in (control) with very low dispersion
#' @return named list of true abundances and observed abundances
#' @export
resample_counts <- function(theta, groups, size_factor_correlation, spike_in = FALSE) {
  p <- (length(theta) - 2)/2
  n_total <- length(groups)
  abundances <- sapply(1:p, function(j) {
    offset <- 2 + (j-1)*2 + 1
    mu.ij <- exp(theta[offset] + groups*theta[offset+1])
    if(spike_in & (j == p)) {
      # rnbinom(n_total, mu = mu.ij, size = 10) # less dispersion
      exp(mu.ij) # noiseless spike-in
    } else {
      rnbinom(n_total, mu = mu.ij, size = theta[1])
    }
  })
  # the abundances matrix has dimensions n*2 samples x p genes
  
  # transform these to proportions
  sampled_proportions <- apply(abundances, 1, function(x) x/sum(x))
  realized_total_counts <- rowSums(abundances)
  size_factors.observed_counts <- round(simcor(realized_total_counts, mean(realized_total_counts), sd(realized_total_counts),
                                               size_factor_correlation))
  # since we're sampling from a normal around the observed counts occasionally we'll get a suggested
  # total count that is negative; reflect this
  # this increases correlation but barely
  size_factors.observed_counts <- abs(size_factors.observed_counts)
  
  observed_counts <- sapply(1:n_total, function(i) {
    if(is.na(size_factors.observed_counts[i])) {
      print(i)
    }
    rmultinom(1, size = size_factors.observed_counts[i], prob = sampled_proportions[,i])
  })
  observed_counts <- t(observed_counts)
  # the observed_counts matrix has dimensions n*2 samples x p genes
  
  return(list(abundances = abundances, observed_counts = observed_counts))
}

#' Simulate UMI count single-cell RNA-seq data
#'
#' @param p number of genes
#' @param n number of samples per group
#' @param k number of distinct cell types
#' @param proportion_da proportion of differentially abundant genes
#' @param size_factor_correlation correlation of observed abundances to original total abundances
#' @param spike_in simulate a control spike in with low dispersion
#' @param possible_fold_changes if specified, a distribution of possible fold changes available for differential abundance simulation
#' @return named list of true abundances, observed abundances, and group assignments
#' @details Samples the mean gene abundances from Halpern et al. (2017) as a references. This gives a data set
#' with on average 75% zeros which may still not be representative in terms of sparsity. Need to follow up. (7/22/2020)
#' @import LaplacesDemon
#' @export
simulate_singlecell_RNAseq <- function(p = 20000, n = 500, k = 1, proportion_da = 0.1, size_factor_correlation = 0, spike_in = FALSE, possible_fold_changes = NULL) {
  proportions_components <- rdirichlet(1, alpha = rep(1/k, k)*5)
  counts_components <- round(proportions_components[1:(k-1)]*n)
  counts_components <- c(counts_components, n - sum(counts_components))
  # if differentially expressed, a gene will be differentially expressed in all cell types (for now)
  Halpern2017_data <- readRDS("data/summary_Halpern2017_base.rds")
  sampled_relative_abundances <- sample(Halpern2017_data$mean_relative_abundances, size = p, replace = TRUE)
  # re-close
  sampled_relative_abundances <- sampled_relative_abundances / sum(sampled_relative_abundances)
  # convert these to baseline abundances
  baseline_abundances <- sampled_relative_abundances * 20000
  # this grossly looks like real data
  #   plot(density(log(round(baseline_abundances) + 0.5)))
  # assign control vs. treatment groups
  groups <- c(rep(0, n), rep(1, n))
  
  n_da <- round(p * proportion_da)
  da_genes <- as.logical(rbinom(p, size = 1, prob = proportion_da))
  # pick a random change for each differential gene
  da_factor <- rep(1, p)
  
  if(!is.null(possible_fold_changes)) {
    # method 1: randomly selection some proportion of genes to be differentially abundant between conditions
    da_factor[da_genes] <- sample(c(1/possible_fold_changes, possible_fold_changes), size = sum(da_genes), replace = TRUE)
  } else {
    # method 2: use empirical data (tissue vs. tissue expression levels) to choose a realistic fold change
    #           given a gene's quantile of expression; in particular, some lowly expressed genes are
    #           observed jumping between modes (on/off) of expression
    GTEx_data <- readRDS("data/GTEx_empirical_DA.rds")
    da_genes_idx <- which(da_genes == TRUE)
    da_factor[da_genes_idx] <- sapply(da_genes_idx, function(x) {
      x_quantile <- sum(baseline_abundances <= baseline_abundances[x])/p
      list_idx <- 1
      if(x_quantile < 0.6) {
        list_idx <- 2
      } else if(x_quantile < 0.7) {
        list_idx <- 3
      } else if(x_quantile < 0.8) {
        list_idx <- 4
      } else if(x_quantile < 0.9) {
        list_idx <- 5
      } else {
        list_idx <- 6
      }
      # Sample from fold changes that are at least +/- a 2-fold change
      # Anything less than this definitely will not be "visible" as differential expression
      fold_change_list <- GTEx_data[[list_idx]]
      fold_change_list <- fold_change_list[fold_change_list >= 5 | fold_change_list <= (1/5)]
      sample(fold_change_list)[1]
    })
  }
  
  # sample true counts (format: samples [rows] x genes [columns])
  abundances <- sapply(1:p, function(j) {
    # simulate different baselines for different cell types (TBD; see below)
    counts <- rnbinom(n*2, mu = c(rep(baseline_abundances[j], n), rep(baseline_abundances[j]*da_factor[j], n)), size = 1)
    # spike-in also TBD
    counts
  })

  realized_total_counts <- rowSums(abundances)
  size_factors.observed_counts <- round(simcor(realized_total_counts, mean(realized_total_counts), sd(realized_total_counts),
                                               size_factor_correlation))
  # since we're sampling from a normal around the observed counts occasionally we'll get a suggested
  # total count that is negative; reflect this
  # this increases correlation but barely
  size_factors.observed_counts <- abs(size_factors.observed_counts)
  
  # transform to proportions
  sampled_proportions <- apply(abundances, 1, function(x) x/sum(x))
  observed_counts <- sapply(1:(n*2), function(i) {
    resample_probs <- sampled_proportions[,i]
    # censor very low concentration stuff to simulate drop-outs
    # I'm trying to simulate some minimum detectable concentration
    # this threshold is pretty arbitrary; we need to tune this to resemble real data
    resample_probs[resample_probs < 1e-8] <- 0
    rmultinom(1, size = size_factors.observed_counts[i], prob = resample_probs)
  })
  observed_counts <- t(observed_counts)
  
  return(list(abundances = abundances,
              observed_counts = observed_counts,
              groups = groups,
              da_genes = which(da_genes == TRUE)))
}

#' Fit negative binomial model to null and full models and evaluate differential abundance for a focal gene
#'
#' @param data simulated data set
#' @param feature_idx index of gene to test for differential abundance
#' @param call_abundances if TRUE, call DA on abundances; if FALSE, on observed counts
#' @return p-value from NB GLM fit with MASS::glm.nb associated with group coefficient
#' @import MASS
#' @export
call_DA_NB <- function(data, feature_idx, call_abundances = TRUE) {
  if(call_abundances) {
    gene_data <- data.frame(counts = data$abundances[,feature_idx], groups = data$groups)
  } else {
    gene_data <- data.frame(counts = data$observed_counts[,feature_idx], groups = data$groups)
  }
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
    # model fit failed
    data_type <- "abundances"
    if(!call_abundances) {
      data_type <- "observed counts"
    }
    cat(paste0("Fit failed on ",data_type," #",feature_idx,"\n"))
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
#' @return p-value associated with group coefficient
#' @import edgeR
#' @export
call_DA_edgeR <- function(data, call_abundances = TRUE) {
  # DGEList expects samples as columns
  if(call_abundances) {
    dge_obj <- DGEList(counts = t(data$abundances), group = factor(data$groups))
  } else {
    dge_obj <- DGEList(counts = t(data$observed_counts), group = factor(data$groups))
  }
  dge_obj <- calcNormFactors(dge_obj, method = "none") # library size normalization
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
#' @param method differential abundance testing method to use; options are "NB", "GLM", "edgeR"
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
        pval.abundances <- call_DA_CODA(logratios.abundances, i, data$groups)
        pval.observed_counts <- call_DA_CODA(logratios.observed_counts, i, data$groups)
        if(!is.na(pval.abundances) & !is.na(pval.observed_counts)) {
          calls.abundances <- c(calls.abundances, pval.abundances <= alpha/length(evaluate_features))
          calls.observed_counts <- c(calls.observed_counts, pval.observed_counts <= alpha/length(evaluate_features))
        }
      }
    }
  } else {
    if(method == "edgeR") {
      # Note: This inherently looks for DE over relative abundances!
      # pval.abundances <- call_DA_edgeR(data, call_abundances = TRUE)
      # calls.abundances <- pval.abundances <= alpha/length(evaluate_features)
      pval.observed_counts <- call_DA_edgeR(data, call_abundances = FALSE)
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
          if(!is.na(pval.abundances) & !is.na(pval.observed_counts)) {
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
#' @param size_factor_correlation correlation of observed abundances to original total abundances
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
run_RNAseq_evaluation_instance <- function(p, n, proportion_da, k = 1, size_factor_correlation = 0,
                                           alpha = 0.05, use_ALR = FALSE, filter_abundance = 1,
                                           methods = c("NB", "edgeR")) {
  data <- simulate_singlecell_RNAseq(p = p, n = n, k = k, proportion_da = proportion_da, size_factor_correlation = size_factor_correlation,
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
    data_obj$size_factor_correlation <- size_factor_correlation
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
                                     size_factor_correlation = sf_corr,
                                     alpha = 0.05, use_ALR = use_ALR, filter_abundance = filter_abundance,
                                     methods = methods)
    }
  }
}

