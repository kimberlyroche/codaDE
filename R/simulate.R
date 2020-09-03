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

#' Simulate differential abundance parameters (NOT IN USE)
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
  counts <- resample_counts(parameters$theta, groups, size_factor_correlation, spike_in)
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
#' @details Samples the mean gene abundances from Halpern et al. (2017) as a references. This gives a data set
#' with on average 75% zeros which may still not be representative in terms of sparsity. Need to follow up. (7/22/2020)
#' @import LaplacesDemon
#' @export
simulate_singlecell_RNAseq <- function(p = 20000, n = 500, k = 1, proportion_da = 0.1, size_factor_correlation = 0, spike_in = FALSE) {
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
  
  # method 1: randomly selection some proportion of genes to be differentially abundant between conditions
  # fold_changes <- c(3, 4, 5)
  # da_factor[da_genes] <- sample(c(1/fold_changes, fold_changes), size = sum(da_genes), replace = TRUE)
  
  # method 2: use empirical data (tissue vs. tissue expression levels) to choose a realistic fold change
  #           given a gene's quantile of expression; in particular, some lowly expressed genes are
  #           observed jumping between modes (on/off) of expression
  GTEx_data <- readRDS("data/GTEx_empirical_DA.rds")
  da_genes_idx <- which(da_genes == TRUE)
  da_factor[da_genes_idx] <- sapply(da_genes_idx, function(x) {
    x_quantile <- sum(baseline_abundances <= baseline_abundances[x])/p
    if(x_quantile < 0.5) {
      sample(GTEx_data[[1]])[1]
    } else if(x_quantile < 0.6) {
      sample(GTEx_data[[2]])[1]
    } else if(x_quantile < 0.6) {
      sample(GTEx_data[[3]])[1]
    } else if(x_quantile < 0.6) {
      sample(GTEx_data[[4]])[1]
    } else if(x_quantile < 0.6) {
      sample(GTEx_data[[5]])[1]
    } else {
      sample(GTEx_data[[6]])[1]
    }
  })
  
  # sample true counts (format: samples [rows] x genes [columns])
  abundances <- sapply(1:p, function(j) {
    # simulate different baselines for different cell types (TBD; see below)
    counts <- rnbinom(n*2, mu = c(rep(baseline_abundances[j], n), rep(baseline_abundances[j]*da_factor[j], n)), size = 1)
    # spike-in also TBD
    counts
  })
  # this grossly looks like real data in terms of log count distribution
  #   plot(density(log(round(abundances) + 0.5)))
  # and again if we look at the profiles of individual genes
  #   counts <- as.matrix(readRDS("data/counts_Halpern2017_base.rds")$counts) # genes (rows) x samples (columns)
  #   emp_log_mean <- rowMeans(log(counts + 0.5))
  #   highly_expressed_idx <- unname(which(emp_log_mean > 3))
  #   lowly_expressed_idx <- unname(which(emp_log_mean < 1))
  #   high_sample_vec <- unname(counts[sample(highly_expressed_idx)[1],])
  #   low_sample_vec <- unname(counts[sample(lowly_expressed_idx)[1],])
  #   par(mfrow = c(2,2))
  #   plot(1:length(high_sample_vec), high_sample_vec)
  #   plot(rnbinom(1000, mu = mean(high_sample_vec), size = 1))
  #   plot(1:length(high_sample_vec), low_sample_vec)
  #   plot(rnbinom(1000, mu = mean(low_sample_vec), size = 1))
  
  # and if we pick a randomly differentially abundant gene this should be somewhat visible
  #   idx <- sample(which(da_genes == TRUE))[1]
  #   da_factor[idx]
  #   plot(abundances[,idx])
  
  # parameters <- build_theta(beta.0, proportion_da, spike_in)
  # abundances <- sapply(1:p, function(j) {
  #   offset <- 2 + (j-1)*2 + 1
  #   mean_log_expr_celltypes <- c()
  #   for(kk in 1:k) {
  #     mean_log_expr_celltypes <- c(mean_log_expr_celltypes, rep(rnorm(1, parameters$theta[offset], 2), counts_components[kk]))
  #   }
  #   counts <- rnbinom(n*2, mu = exp(mean_log_expr_celltypes + groups*parameters$theta[offset+1]), size = parameters$theta[1])
  #   if(spike_in & (j == p)) {
  #     # for now, simulate the same noiseless spike-in
  #     counts <- c(counts, exp(parameters$theta[offset] + groups*parameters$theta[offset+1]))
  #   }
  #   counts
  # })
  
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
  
  # cat(paste0("Percent zeros in data set: ",round(sum(observed_counts == 0)/(nrow(observed_counts)*ncol(observed_counts)), 2)*100,"\n"))
  
  # if we wanted to calculate entropy for a whole data set, here's how we might go about it
    
  # log_means <- log(colMeans(observed_counts[,parameters$da_genes]) + 1)
  # breaks <- seq(from = min(log_means), to = max(log_means), length.out = 10)
  # chunked_data <- cut(log_means, breaks, 1:9)
  # ent_val <- entropy(table(chunked_data))
  
  # the observed_counts matrix has dimensions n*2 samples x p genes
  # return(list(abundances = abundances,
  #             observed_counts = observed_counts,
  #             groups = groups,
  #             theta = parameters$theta,
  #             da_genes = parameters$da_genes))
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

#' #' Write simulation result to a file
#' #'
#' #' @param out_str string to write to file
#' #' @param out_file log file path and name
#' #' @return NULL
#' #' @import filelock
#' #' @export
#' record_result <- function(out_str, out_file) {
#'   lock_file <- paste0(out_file, ".lck")
#'   lck <- lock(lock_file, timeout = Inf)
#'   write(out_str, file = out_file, append = TRUE)
#'   unlock(lck)
#' }

#' Update the metadata file in the simulated data set directory
#'
#' @param filename autogenerated filename for simulated data set
#' @param p number of genes
#' @param n number of samples per group
#' @param proportion_da proportion of differentially abundant genes
#' @param size_factor_correlation correlation of observed abundances to original total abundances
#' @param k number of cell types to simulate; if not NULL, single-cell data is simulated; if NULL, bulk RNA-seq data is simulated
#' @return NULL
#' @import filelock
#' @export
update_metadata <- function(filename, p, n, proportion_da, size_factor_correlation, k) {
  append_row <- data.frame(filename = filename,
                           p = p,
                           n = n,
                           proportion_da = proportion_da,
                           size_factor_correlation = size_factor_correlation,
                           k = k)
  metadata_file <- file.path("simulated_data", "metadata.tsv") # hard-code this file path for now
  if(!file.exists(metadata_file)) {
    # create it
    write.table(append_row, file = metadata_file, quote = FALSE, sep = '\t', row.names = FALSE)
  } else {
    lock_file <- paste0(metadata_file, ".lck")
    lck <- lock(lock_file, timeout = Inf)
    metadata <- read.table(metadata_file, header = T, stringsAsFactors = FALSE)
    metadata <- rbind(metadata, append_row)
    write.table(metadata, file = metadata_file, quote = FALSE, sep = '\t', row.names = FALSE)
    dump <- unlock(lck)
  }
}

#' Generate simulated data with RNA-seq like features and > 3 fold differential abundance, evaluate error in
#' differential abundance calls and write these to an output file
#'
#' @param data simulated data set
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
evaluate_DA <- function(data, alpha = 0.05, use_ALR = FALSE, filter_abundance = 1, call_DA_by_NB = TRUE, rarefy = FALSE) {
  # calculate and compare differential abundance calls
  cat("Evaluating differential abundance...\n")
  calls.abundances <- c()
  calls.observed_counts <- c()

  if(rarefy) {
    rarefy_total <- min(apply(data$observed_counts, 1, sum))
  } else {
    rarefy_total <- 0
  }  
  evaluate_features <- apply(data$observed_counts, 2, function(x) mean(x) > filter_abundance)

  p <- ncol(data$abundances)
  
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
          calls.abundances <- c(calls.abundances, pval.abundances <= alpha/length(evaluate_features))
          calls.observed_counts <- c(calls.observed_counts, pval.observed_counts <= alpha/length(evaluate_features))
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

  # calculate the log mean expression for the DA group at baseline and include a measure of how different
  #     this is from overall log mean expression at baseline -- i.e. is our selection of DE genes biased?
  baseline_expr_all <- colMeans(log(data$abundances[data$groups == 0,] + 0.5))
  detected_da <- which(calls.abundances == TRUE)
  if(length(detected_da) > 0) {
    baseline_expr_da <- colMeans(log(data$abundances[data$groups == 0, detected_da] + 0.5))
    ordered_baseline <- baseline_expr_all[order(baseline_expr_all)]
    median_expr_da_quantile <- sum(ordered_baseline < median(baseline_expr_da))/length(ordered_baseline)

    # calculate the log mean expression for the DA group under "treatment"; the median of the signed different
    #     will be a predictor too
    affected_expr_dat <- colMeans(log(data$abundances[data$groups == 1, detected_da] + 0.5))
    net_dir_da <- median(affected_expr_dat - baseline_expr_da)
  } else {
    median_expr_da_quantile <- NULL
    net_dir_da <- NULL
  }

  # sparsity level (percent zeros)
  sparsity <- sum(data$observed_counts == 0)/(nrow(data$observed_counts)*ncol(data$observed_counts))
  
  # uniformity of composition at baseline
  binned_log_expression <- as.numeric(table(cut(baseline_expr_all, breaks = seq(from = min(baseline_expr_all), to = max(baseline_expr_all), length.out = 10))))
  max_entropy <- entropy(rep(round(sum(binned_log_expression) / length(binned_log_expression)), length(binned_log_expression)))
  sim_entropy <- entropy(binned_log_expression) / max_entropy

  save_obj <- list(data = data,
                   properties = list(median_expr_da_quantile = median_expr_da_quantile,
                                     net_dir_da = net_dir_da,
                                     sparsity = sparsity,
                                     sim_entropy = sim_entropy
                   ),
                   results = list(FN = FN,
                                  FP = FP,
                                  TN = TN,
                                  TP = TP,
                                  TPR = TP / (TP + FN),
                                  FPR = FP / (FP + TN),
                                  quantity_evaluated = quantity_evaluated,
                                  call_DA_by_NB = call_DA_by_NB,
                                  filter_abundance = filter_abundance,
                                  no_features_above_threshold = sum(evaluate_features),
                                  no_features_detectable = sum(calls.abundances == TRUE),
                                  calls.abundances = calls.abundances,
                                  calls.observed_counts = calls.observed_counts))

  return(save_obj)
}

#' Generate simulated data with RNA-seq like features and > 3 fold differential abundance, evaluate error in
#' differential abundance calls and write these to an output file
#'
#' @param p number of genes
#' @param n number of samples per group
#' @param proportion_da proportion of differentially abundant genes
#' @param k number of cell types to simulate; if not NULL, single-cell data is simulated; if NULL, bulk RNA-seq data is simulated
#' @param size_factor_correlation correlation of observed abundances to original total abundances
#' @param alpha significant level below which to call a feature differentially abundant
#' @param use_ALR if TRUE, evaluates differential abundance using a spike-in and the additive logratio
#' @param filter_abundance the minimum average abundance to require of features we'll evaluate for differential abundance
#' (e.g. a filter_abundance of 1 evaluates features with average abundance across conditions of at least 1)
#' @param call_DA_by_NB if TRUE, uses a NB GLM to call differential abundance; if FALSE, uses a log linear model with a
#' permutation test
#' @param rarefy if TRUE, resamples the counts in each sample to the lowest observed total counts
#' @param analysis_label this is a short name for this analysis that should correspond to a folder name in simulated_analyses
#' @details Writes out error and simulation run statistics to a file.
#' @return NULL
#' @import entropy
#' @import uuid
#' @export
run_RNAseq_evaluation_instance <- function(p, n, proportion_da, k = NULL, size_factor_correlation = 0,
                                           alpha = 0.05, use_ALR = FALSE, filter_abundance = 1,
                                           call_DA_by_NB = TRUE, rarefy = FALSE, analysis_label = NULL) {
  if(is.null(analysis_label)) {
    stop("Missing analysis label!\n")
  }
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
  
  save_obj <- evaluate_DA(data, alpha = alpha, use_ALR = use_ALR, filter_abundance = filter_abundance, call_DA_by_NB = call_DA_by_NB, rarefy = rarefy)
  data_obj <- save_obj
  data_obj$properties$p <- p
  data_obj$properties$proportion_da <- proportion_da
  data_obj$properties$size_factor_correlation <- size_factor_correlation
  data_obj$results <- NULL
  results_obj <- save_obj$results
  new_filename <- paste0(UUIDgenerate(), ".rds")
  saveRDS(data_obj, file.path("simulated_data", new_filename))
  saveRDS(results_obj, file.path("simulated_analyses", analysis_label, new_filename))
  update_metadata(new_filename, p, n, proportion_da, size_factor_correlation, k)
}

#' Generate simulated data with RNA-seq like features and > 3 fold differential abundance, evaluate error in
#' differential abundance calls and write these to an output file
#'
#' @param p number of genes
#' @param n number of samples per group
#' @param proportion_da proportion of differentially abundant genes
#' @param k number of cell types to simulate; if not NULL, single-cell data is simulated; if NULL, bulk RNA-seq data is simulated
#' @param size_factor_correlation correlation of observed abundances to original total abundances
#' @param alpha significant level below which to call a feature differentially abundant
#' @param use_ALR if TRUE, evaluates differential abundance using a spike-in and the additive logratio
#' @param filter_abundance the minimum average abundance to require of features we'll evaluate for differential abundance
#' (e.g. a filter_abundance of 1 evaluates features with average abundance across conditions of at least 1)
#' @param call_DA_by_NB if TRUE, uses a NB GLM to call differential abundance; if FALSE, uses a log linear model with a
#' permutation test
#' @param rarefy if TRUE, resamples the counts in each sample to the lowest observed total counts
#' @param analysis_label this is a short name for this analysis that should correspond to a folder name in simulated_analyses
#' @details Writes out error and simulation run statistics to a file.
#' @return NULL
#' @import filelock
#' @export
evaluate_existing_RNAseq_instance <- function(p, proportion_da, k, size_factor_correlation,
                                              alpha = 0.05, use_ALR = FALSE, filter_abundance = 1,
                                              call_DA_by_NB = TRUE, rarefy = FALSE, analysis_label = NULL) {
  metadata_file <- file.path("simulated_data", "metadata.tsv") # assume this isn't being actively written to!
  metadata <- read.table(metadata_file, header = T, stringsAsFactors = FALSE)
  # find filenames matching conditions
  usable_filenames <- metadata[metadata$p == p &
                               metadata$proportion_da == proportion_da &
                               metadata$k == k &
                               metadata$size_factor_correlation == size_factor_correlation,]$filename

  # check the ledger to see which are not in use
  ledger_file <- file.path("simulated_analyses", analysis_label, "ledger.tsv")
  if(!file.exists(ledger_file)) {
    # create it
    file_to_use <- sample(usable_filenames)[1]
    write.table(data.frame(filename = file_to_use), file = ledger_file, quote = FALSE, sep = '\t', row.names = FALSE)
  } else {
    lock_file <- paste0(ledger_file, ".lck")
    lck <- lock(lock_file, timeout = Inf)
    ledger <- read.table(ledger_file, header = T, stringsAsFactors = FALSE)
    file_to_use <- sample(usable_filenames[!(usable_filenames %in% ledger$filename)])[1]
    if(!is.na(file_to_use)) {
      ledger <- rbind(ledger, data.frame(filename = file_to_use))
      write.table(ledger, file = ledger_file, quote = FALSE, sep = '\t', row.names = FALSE)
    }
    dump <- unlock(lck)
  }

  if(!is.na(file_to_use)) {
    # run the differential abundance test
    cat("Evaluating file",file_to_use,"\n")
    data_filename <- file.path("simulated_data", file_to_use)
    data <- readRDS(data_filename)
    save_obj <- evaluate_DA(data$data, alpha = alpha, use_ALR = use_ALR, filter_abundance = filter_abundance, call_DA_by_NB = call_DA_by_NB, rarefy = rarefy)
    saveRDS(save_obj$results, file.path("simulated_analyses", analysis_label, file_to_use))
  }
}

#' Generate simulated data with RNA-seq like features and > 3 fold differential abundance, evaluate false negatives and positives in
#' differential abundance calls and write these to an output file
#'
#' @param p number of genes
#' @param n number of samples per group
#' @param k number of cell types to simulate; if not NULL, single-cell data is simulated; if NULL, bulk RNA-seq data is simulated
#' @param de_sweep vector of proportions of differentially abundant genes to simulate
#' @param corr_sweep vector of correlations of size factors (true vs. observed abundances) to simulate
#' @param alpha significant level below which to call a feature differentially abundant
#' @param use_ALR if TRUE, evaluates differential abundance using a spike-in and the additive logratio
#' @param filter_abundance the minimum average abundance to require of features we'll evaluate for differential abundance
#' (e.g. a filter_abundance of 1 evaluates features with average abundance across conditions of at least 1)
#' @param call_DA_by_NB if TRUE, uses a NB GLM to call differential abundance; if FALSE, uses a log linear model with a
#' permutation test
#' @param rarefy if TRUE, resamples the counts in each sample to the lowest observed total counts
#' @param use_existing_simulations if TRUE looks for existing simulations (specified in a metadata file) rather than simulating new data sets
#' @param analysis_label this is a short name for this analysis that should correspond to a folder name in simulated_analyses
#' @details Writes out error and simulation run statistics to a file.
#' @return NULL
#' @export
sweep_simulations <- function(p, n, k = NULL, de_sweep = seq(from = 0.1, to = 0.9, by = 0.1),
                         corr_sweep = seq(from = 0.1, to = 0.9, by = 0.1), alpha = 0.05, use_ALR = FALSE,
                         filter_abundance = 0, call_DA_by_NB = TRUE, rarefy = FALSE,
                         use_existing_simulations = FALSE, analysis_label = NULL) {
  if(is.null(analysis_label)) {
    stop("Missing analysis label!\n")
  }
  for(de_prop in de_sweep) {
    for(sf_corr in corr_sweep) {
      out_str <- paste0("w/ ",p," genes, DA proportion = ",round(de_prop, 2),", and size factor correlation = ",round(sf_corr, 2),"\n")
      if(is.null(k)) {
        cat("Bulk RNA-seq",out_str)
      } else {
        cat("Single-cell RNA-seq",out_str)
      }
      if(use_existing_simulations) {
        evaluate_existing_RNAseq_instance(p, proportion_da = de_prop, k = k, size_factor_correlation = sf_corr,
                                                      alpha = 0.05, use_ALR = use_ALR, filter_abundance = filter_abundance,
                                                      call_DA_by_NB = call_DA_by_NB, rarefy = rarefy, analysis_label = analysis_label)
      } else {
        run_RNAseq_evaluation_instance(p, n, proportion_da = de_prop, k = k,
                                       size_factor_correlation = sf_corr,
                                       alpha = 0.05, use_ALR = use_ALR, filter_abundance = filter_abundance,
                                       call_DA_by_NB = call_DA_by_NB, rarefy = rarefy, analysis_label = analysis_label)
      }
    }
  }
}









