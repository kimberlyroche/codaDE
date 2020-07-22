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
#' @return named list of true abundances and observed abundances
#' @export
resample_counts <- function(theta, groups) {
  p <- (length(theta) - 2)/2
  abundances <- sapply(1:p, function(j) {
    offset <- 2 + (j-1)*2 + 1
    mu.ij <- exp(theta[offset] + groups*theta[offset+1])
    if(spike_in & (j == p)) {
      # rnbinom((n*2), mu = mu.ij, size = 10) # less dispersion
      exp(mu.ij) # noiseless spike-in
    } else {
      rnbinom((n*2), mu = mu.ij, size = theta[1])
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
  
  observed_counts <- sapply(1:(n*2), function(i) {
    if(is.na(size_factors.observed_counts[i])) {
      print(i)
    }
    rmultinom(1, size = size_factors.observed_counts[i], prob = sampled_proportions[,i])
  })
  observed_counts <- t(observed_counts)
  # the observed_counts matrix has dimensions n*2 samples x p genes
  
  return(list(abundances = abundances, observed_counts = observed_counts))
}

#' Simulate differential abundance parameters
#'
#' @param beta.0 vector of log baseline (mean) expression for all genes
#' @param proportion_da proportion of differentially abundant genes
#' @param spike_in flag indicating whether or not to simulate a spike-in (control) with very low dispersion
#' @return complete vector of dispersion and log coefficient parameters
#' @export
build_theta <- function(beta.0, proportion_da, spike_in) {
  p <- length(beta.0)
  theta <- rep(0.25, 2) # fix the dispersions between groups for simplicity now but we should consider
  # varying this as well
  if(spike_in) {
    # simulate a spike-in at the median log abundance and tack it on
    spike_in_beta <- median(beta.0)
    beta.0 <- c(beta.0, spike_in_beta)
  }
  # select a random set of features to perturb
  da_genes <- sort(sample(1:p)[1:round(proportion_da*p)])
  # simulate a random large fold change (between +/- 3 and 5) for these, otherwise give the group effect (beta.1)
  # as zero
  beta.1 <- sapply(1:p, function(idx) {
    if(idx %in% da_genes) {
      fc <- sample(c(3:5))[1]
      fc_sign <- 1
      if(rbinom(1, size = 1, prob = 0.5)) {
        fc_sign <- -1
      }
      fc_sign*log(fc)
    } else {
      0
    }
  })
  # the set of (2) group dispersion parameters and (2) per-gene log covariates is theta
  theta <- c(theta, c(rbind(beta.0, beta.1)))
  return(list(theta = theta, da_genes = da_genes))
}

#' Simulate bulk RNA-seq-like data
#'
#' @param p number of genes
#' @param n number of samples per group (e.g. control vs. treatment)
#' @param proportion_da proportion of differentially abundant genes
#' @param size_factor_correlation correlation of observed abundances to original absolute abundances
#' @param spike_in flag indicating whether or not to simulate a spike-in (control) with very low dispersion
#' @return named list of true abundances, observed abundances, and group assignments
#' @export
simulate_bulk_RNAseq <- function(p = 20000, n = 500, proportion_da = 0.1, size_factor_correlation = 0, spike_in = FALSE) {
  data_dir <- "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"
  GTEx_stats <- readRDS(file.path(data_dir, "empirical_sample.rds"))
  # sample a set of baseline log abundances for p genes
  beta.0 <- sample(log(GTEx_stats$mean_abundance_profiles + 0.5), size = p, replace = TRUE)
  groups <- c(rep(0, n), rep(1, n))
  parameters <- build_theta(beta.0, proportion_da, spike_in)
  counts <- resample_counts(parameters$theta, groups)
  return(list(abundances = counts$abundances,
              observed_counts = counts$observed_counts,
              groups = groups,
              theta = parameters$theta,
              da_genes = parameters$da_genes))
}

#' Simulate UMI count single-cell RNA-seq data
#'
#' @param p number of genes
#' @param n number of samples per group
#' @param k number of distinct cell types
#' @param proportion_da proportion of differentially abundant genes
#' @param size_factor_correlation correlation of observed abundances to original total abundances
#' @param spike_in simulate a control spike in with low dispersion
#' @return named list of true abundances, observed abundances, and group assignments
#' @import LaplacesDemon
#' @export
simulate_singlecell_RNAseq <- function(p = 20000, n = 500, k = 1, proportion_da = 0.1, size_factor_correlation = 0, spike_in = FALSE) {
  proportions_components <- rdirichlet(1, alpha = rep(1/k, k)*5)
  counts_components <- round(proportions_components[1:(k-1)]*n)
  counts_components <- c(counts_components, n - sum(counts_components))
  # if differentially expressed, a gene will be differentially expressed in all cell types (for now)
  beta.0 <- rnorm(p, log(30), 1)
  groups <- c(rep(0, n), rep(1, n))
  parameters <- build_theta(beta.0, proportion_da, spike_in)
  abundances <- sapply(1:p, function(j) {
    offset <- 2 + (j-1)*2 + 1
    mean_log_expr_celltypes <- c()
    for(kk in 1:k) {
      mean_log_expr_celltypes <- c(mean_log_expr_celltypes, rep(rnorm(1, parameters$theta[offset], 2), counts_components[kk]))
    }
    counts <- rnbinom(n*2, mu = exp(mean_log_expr_celltypes + groups*parameters$theta[offset+1]), size = parameters$theta[1])
    if(spike_in & (j == p)) {
      # for now, simulate the same noiseless spike-in
      counts <- c(counts, exp(parameters$theta[offset] + groups*parameters$theta[offset+1]))
    }
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
  # the observed_counts matrix has dimensions n*2 samples x p genes
  return(list(abundances = abundances,
              observed_counts = observed_counts,
              groups = groups,
              theta = parameters$theta,
              da_genes = parameters$da_genes))
}

#' Fit negative binomial model to null and full models and evaluate differential abundance for a focal gene
#'
#' @param data simulated data set
#' @param feature_idx index of gene to test for differential abundance
#' @param call_abundances if TRUE, call DA on abundances; if FALSE, on observed counts
#' @param rarefy_depth if non-zero, specifies the depth at which to resample the counts
#' @return p-value from NB GLM fit with MASS::glm.nb
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
    return(gene_data)
  }
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
#' @param rarefy if TRUE, resamples the counts in each sample to the lowest observed total counts
#' @details Writes out error and simulation run statistics to a file.
#' @return NULL
#' @export
run_RNAseq_evaluation_instance <- function(p, n, proportion_da, run_label, k = NULL, size_factor_correlation = 0,
                                           output_file = "results.txt", alpha = 0.05, use_ALR = FALSE,
                                           filter_abundance = 1, rarefy = FALSE) {
  all_taxa_observed <- FALSE
  while(!all_taxa_observed) {
    cat("Simulating data...\n")
    if(is.null(k)) {
      data <- simulate_bulk_RNAseq(p = p, n = n, proportion_da = proportion_da, size_factor_correlation = size_factor_correlation,
                                   spike_in = use_ALR)
    } else {
      data <- simulate_singlecell_RNAseq(p = p, n = n, k = k, proportion_da = proportion_da, size_factor_correlation = size_factor_correlation,
                                         spike_in = use_ALR)
    }
    g0.ab <- data$abundances[data$groups == 0,] # samples x genes
    g1.ab <- data$abundances[data$groups == 1,] # samples x genes
    g0.oc <- data$observed_counts[data$groups == 0,] # samples x genes
    g1.oc <- data$observed_counts[data$groups == 1,] # samples x genes
    if(sum(rowSums(g0.ab) == 0) == 0 & sum(colSums(g0.ab) == 0) == 0 &
       sum(rowSums(g1.ab) == 0) == 0 & sum(colSums(g1.ab) == 0) == 0 &
       sum(rowSums(g0.oc) == 0) == 0 & sum(colSums(g0.oc) == 0) == 0 &
       sum(rowSums(g1.oc) == 0) == 0 & sum(colSums(g1.oc) == 0) == 0) {
      all_taxa_observed <- TRUE
    }
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
      if(i %% 1000 == 0) {
        cat("Evaluating DA on feature:",i,"\n")
      }
      if(evaluate_features[i]) {
        pval.abundances <- call_DA(data, i, call_abundances = TRUE, rarefy = rarefy_total)
        pval.observed_counts <- call_DA(data, i, call_abundances = FALSE, rarefy = rarefy_total)
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
  out_str <- paste0(run_label,"\t",
                    p,"\t",
                    proportion_da,"\t",
                    size_factor_correlation,"\t",
                    quantity_evaluated,"\t",
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
#' @param rarefy if TRUE, resamples the counts in each sample to the lowest observed total counts
#' @details Writes out error and simulation run statistics to a file.
#' @return NULL
#' @export
sweep_simulations <- function(p, n, run_label, k = NULL, de_sweep = seq(from = 0.1, to = 0.9, by = 0.1),
                         corr_sweep = seq(from = 0.1, to = 0.9, by = 0.1), output_file = "results.txt",
                         alpha = 0.05, use_ALR = FALSE, filter_abundance = 0, rarefy = FALSE) {
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
                                     rarefy = rarefy)
    }
  }
}









