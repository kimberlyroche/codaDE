
#' Fit negative binomial model to null and full models and evaluate differential abundance for a focal gene
#'
#' @param data simulated data set
#' @param feature_idx index of gene to test for differential abundance
#' @param call_abundances if TRUE, call DA on abundances; if FALSE, on observed counts
#' @param rarefy_depth if non-zero, specifies the depth at which to resample the counts
#' @return p-value from NB GLM fit with MASS::glm.nb associated with group coefficient
#' @import MASS
#' @export
call_DA <- function(data, feature_idx, call_abundances = TRUE, rarefy_depth = 0) {
  if(call_abundances) {
    gene_data <- data.frame(counts = data$abundances[,feature_idx], groups = data$groups)
  } else {
    gene_data <- data.frame(counts = data$observed_counts[,feature_idx], groups = data$groups)
  }
  if(rarefy_depth > 0) {
    gene_data$counts <- round(gene_data$counts * (rarefy_depth / sum(gene_data$counts)))
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

#' #' Fit log linear model to null and full models and evaluate differential abundance for a focal gene
#' #' This requires a null distribution to have been constructed previously
#' #' 
#' #' @param data simulated data set
#' #' @param use_abundances if TRUE, use abundances to calculate null distro; if FALSE, use observed counts
#' #' @param n_permutations number of times to permute the full data set
#' #' @return p-value associated with group coefficient
#' #' @export
#' generate_null_distribution <- function(data, use_abundances = TRUE, n_permutations = 10) {
#'   null_distribution <- c()
#'   for(p in 1:n_permutations) {
#'     n_samples <- length(data$groups)
#'     groups <- rbinom(n_samples, size = 1, prob = rep(0.5, n_samples))
#'     obs_rsq_vector <- sapply(1:ncol(data$abundances), function(z) {
#'       call_DA_LM_sub(data, z, call_abundances = use_abundances)
#'     })
#'     null_distribution <- c(null_distribution, obs_rsq_vector)
#'   }
#'   null_distribution <- null_distribution[order(null_distribution)]
#'   return(null_distribution)
#' }

#' #' Fit log linear model to null and full models and evaluate differential abundance for a focal gene
#' #' 
#' #' @param data simulated data set
#' #' @param feature_idx index of gene to test for differential abundance
#' #' @param call_abundances if TRUE, call DA on abundances; if FALSE, on observed counts
#' #' @return p-value associated with group coefficient
#' #' @export
#' call_DA_LM_sub <- function(data, feature_idx, call_abundances = TRUE) {
#'   if(call_abundances) {
#'     gene_data <- data.frame(log_counts = log(data$abundances[,feature_idx] + 0.5), groups = data$groups)
#'   } else {
#'     gene_data <- data.frame(log_counts = log(data$observed_counts[,feature_idx] + 0.5), groups = data$groups)
#'   }
#'   fit <- lm(log_counts ~ groups, data = gene_data)
#'   obs_rsq <- summary(fit)$adj.r.squared
#'   return(obs_rsq)
#' }

#' #' Fit log linear model to null and full models and evaluate differential abundance for a focal gene;
#' #' this requires a null distribution to have been constructed previously such that p-values can be
#' #' simulated
#' #' 
#' #' @param data simulated data set
#' #' @param feature_idx index of gene to test for differential abundance
#' #' @param call_abundances if TRUE, call DA on abundances; if FALSE, on observed counts
#' #' @param null_distribution simulated null distribution of R-squared values
#' #' @return p-value associated with group coefficient
#' #' @export
#' call_DA_LM <- function(data, feature_idx, call_abundances = TRUE, null_distribution) {
#'   obs_rsq <- call_DA_LM_sub(data, feature_idx, call_abundances)
#'   pval <- sum(null_distribution > obs_rsq)/length(null_distribution)
#'   return(pval)
#' }

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

#' Write simulation result to a file
#'
#' @param out_str string to write to file
#' @param out_file log file path and name
#' @return NULL
#' @import filelock
#' @export
record_result <- function(out_str, out_file) {
  lock_file <- paste0(out_file, ".lck")
  lck <- lock(lock_file, timeout = Inf)
  write(out_str, file = out_file, append = TRUE)
  unlock(lck)
}

#' Generate simulated data with RNA-seq like features and > 3 fold differential abundance, evaluate error in
#' differential abundance calls and write these to an output file
#'
#' @param p number of genes
#' @param n number of samples per group
#' @param proportion_da proportion of differentially abundant genes
#' @param run_label human-readable string identifier for this simulation run
#' @param k number of cell types to simulate; if not NULL, single-cell data is simulated; if NULL, bulk RNA-seq data is simulated
#' @param size_factor_correlation correlation of observed abundances to original total abundances
#' @param output_file output file to append to
#' @param alpha significant level below which to call a feature differentially abundant
#' @param use_ALR if TRUE, evaluates differential abundance using a spike-in and the additive logratio
#' @param filter_abundance the minimum average abundance to require of features we'll evaluate for differential abundance
#' (e.g. a filter_abundance of 1 evaluates features with average abundance across conditions of at least 1)
#' @param call_DA_by_NB if TRUE, uses a NB GLM to call differential abundance; if FALSE, uses a log linear model with a
#' permutation test
#' @param rarefy if TRUE, resamples the counts in each sample to the lowest observed total counts
#' @details Writes out error and simulation run statistics to a file.
#' @return NULL
#' @export
run_RNAseq_evaluation_instance <- function(p, n, proportion_da, run_label, k = NULL, size_factor_correlation = 0,
                                           output_file = "results.txt", alpha = 0.05, use_ALR = FALSE,
                                           filter_abundance = 1, call_DA_by_NB = TRUE, rarefy = FALSE) {
  cat("Simulating data...\n")
  if(is.null(k)) {
    data <- simulate_bulk_RNAseq(p = p, n = n, proportion_da = proportion_da, size_factor_correlation = size_factor_correlation,
                                 spike_in = use_ALR)
  } else {
    data <- simulate_singlecell_RNAseq(p = p, n = n, k = k, proportion_da = proportion_da, size_factor_correlation = size_factor_correlation,
                                       spike_in = use_ALR)
  }
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
  
  if(rarefy) {
    rarefy_total <- min(apply(data$observed_counts, 1, sum))
  } else {
    rarefy_total <- 0
  }
  
  # calculate and compare differential abundance calls
  cat("Evaluating differential abundance...\n")
  calls.abundances <- c()
  calls.observed_counts <- c()
  
  evaluate_features <- apply(data$observed_counts, 2, function(x) mean(x) > filter_abundance)
  
  if(use_ALR) {
    # this setup is currently only possible with NB DA calling
    logratios.abundances <- get_logratios(data, call_abundances = TRUE)
    logratios.observed_counts <- get_logratios(data, call_abundances = FALSE)
    for(i in 1:p) {
      if(evaluate_features[i]) {
        pval.abundances <- call_DA_CODA(logratios.abundances, i, data$groups)
        pval.observed_counts <- call_DA_CODA(logratios.observed_counts, i, data$groups)
        if(!is.na(pval.abundances) & !is.na(pval.observed_counts)) {
          calls.abundances <- c(calls.abundances, pval.abundances <= alpha/p)
          calls.observed_counts <- c(calls.observed_counts, pval.observed_counts <= alpha/p)
        }
      }
    }
  } else {
    for(i in 1:p) {
      if(i %% 100 == 0) {
        cat("Evaluating DA on feature:",i,"\n")
      }
      if(evaluate_features[i]) {
        if(call_DA_by_NB) {
          pval.abundances <- call_DA(data, i, call_abundances = TRUE, rarefy = rarefy_total)
          pval.observed_counts <- call_DA(data, i, call_abundances = FALSE, rarefy = rarefy_total)
        } else {
          pval.abundances <- call_DA_LM(data, i, call_abundances = TRUE)
          pval.observed_counts <- call_DA_LM(data, i, call_abundances = FALSE)
        }
        if(!is.na(pval.abundances) & !is.na(pval.observed_counts)) {
          calls.abundances <- c(calls.abundances, pval.abundances <= alpha/p)
          calls.observed_counts <- c(calls.observed_counts, pval.observed_counts <= alpha/p)
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
  
  # calculate a couple of quantities that might help us model error in low-feature-number, high DE cases
  # (1) calculate the log mean expression for the DA group at baseline and include a measure of how different
  #     this is from overall log mean expression at baseline -- i.e. is our selection of DE genes biased?
  baseline_expr_all <- colMeans(log(data$abundances[data$groups == 0,] + 0.5))
  baseline_expr_da <- colMeans(log(data$abundances[data$groups == 0, data$da_genes] + 0.5))
  ordered_baseline <- baseline_expr_all[order(baseline_expr_all)]
  median_expr_da_quantile <- sum(ordered_baseline < median(baseline_expr_da))/length(ordered_baseline)
  
  # (2) calculate the log mean expression for the DA group under "treatment"; the median of the signed different
  #     will be a predictor too
  med_affected_expr_dat <- colMeans(log(data$abundances[data$groups == 1, data$da_genes] + 0.5))
  net_dir_da <- med_affected_expr_dat - med_baseline_expr_da
  
  out_str <- paste0(run_label,"\t",
                    p,"\t",
                    proportion_da,"\t",
                    size_factor_correlation,"\t",
                    round(median_expr_da_quantile, 5),"\t",
                    round(median(net_dir_da), 5),"\t",
                    quantity_evaluated,"\t",
                    call_DA_by_NB,"\t",
                    filter_abundance,"\t",
                    TP,"\t",
                    FP,"\t",
                    TN,"\t",
                    FN,"\t",
                    sum(evaluate_features),"\t")
  record_result(out_str, file.path("output",output_file))
}

#' Generate simulated data with RNA-seq like features and > 3 fold differential abundance, evaluate false negatives and positives in
#' differential abundance calls and write these to an output file
#'
#' @param p number of genes
#' @param n number of samples per group
#' @param run_label human-readable string identifier for this simulation run
#' @param k number of cell types to simulate; if not NULL, single-cell data is simulated; if NULL, bulk RNA-seq data is simulated
#' @param de_sweep vector of proportions of differentially abundant genes to simulate
#' @param corr_sweep vector of correlations of size factors (true vs. observed abundances) to simulate
#' @param output_file output file to append to
#' @param alpha significant level below which to call a feature differentially abundant
#' @param use_ALR if TRUE, evaluates differential abundance using a spike-in and the additive logratio
#' @param filter_abundance the minimum average abundance to require of features we'll evaluate for differential abundance
#' (e.g. a filter_abundance of 1 evaluates features with average abundance across conditions of at least 1)
#' @param call_DA_by_NB if TRUE, uses a NB GLM to call differential abundance; if FALSE, uses a log linear model with a
#' permutation test
#' @param rarefy if TRUE, resamples the counts in each sample to the lowest observed total counts
#' @details Writes out error and simulation run statistics to a file.
#' @return NULL
#' @export
sweep_simulations <- function(p, n, run_label, k = NULL, de_sweep = seq(from = 0.1, to = 0.9, by = 0.1),
                              corr_sweep = seq(from = 0.1, to = 0.9, by = 0.1), output_file = "results.txt",
                              alpha = 0.05, use_ALR = FALSE, filter_abundance = 0, call_DA_by_NB = TRUE, rarefy = FALSE) {
  for(de_prop in de_sweep) {
    for(sf_corr in corr_sweep) {
      out_str <- paste0("simulated evaluating size factor correlation = ",round(sf_corr, 2)," and DA proportion = ",round(de_prop, 2),"\n")
      if(is.null(k)) {
        cat("Bulk RNA-seq",out_str)
      } else {
        cat("Single-cell RNA-seq",out_str)
      }
      run_RNAseq_evaluation_instance(p, n, proportion_da = de_prop, run_label = run_label, k = k,
                                     size_factor_correlation = sf_corr, output_file = output_file,
                                     alpha = 0.05, use_ALR = use_ALR, filter_abundance = filter_abundance,
                                     call_DA_by_NB = call_DA_by_NB, rarefy = rarefy)
    }
  }
}
