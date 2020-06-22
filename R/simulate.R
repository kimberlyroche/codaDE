# utility function to simulate a vector with some desired correlation to a reference vector
# source: https://stackoverflow.com/questions/47775389/how-to-simulate-a-vector-that-is-correlated-in-a-different-way-to-two-other-ex
simcor <- function(x, ymean = 0, ysd = 1, correlation = 0) {
  n <- length(x)
  y <- rnorm(n)
  z <- correlation * scale(x)[,1] + sqrt(1 - correlation^2) * scale(resid(lm(y ~ x)))[,1]
  yresult <- ymean + ysd * z
  yresult
}

#' Simulate RNA-seq (like) data
#'
#' @param p number of genes
#' @param n number of samples per group
#' @param proportion_de proportion of differentially expressed genes
#' @param size_factor_correlation correlation of observed abundances to original total abundances
#' @return named list of true abundances, observed abundances, and group assignments
#' @export
simulate_RNAseq <- function(p = 5000, n = 500, proportion_de = 0.5, size_factor_correlation = 0) {
  data_dir <- "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct"
  GTEx_stats <- readRDS(file.path(data_dir, "empirical_sample.rds"))
  groups <- c(rep(0, n), rep(1, n))
  theta <- rep(0.25, 2) # fix the dispersions between groups for simplicity now
  beta.0 <- sample(log(GTEx_stats$mean_abundance_profiles + 0.5), size = p, replace = TRUE)
  de_elements <- sort(sample(1:p)[1:round(proportion_de*p)])
  # simulate a random fold change between (+/-) 3 and 5 for these
  beta.1 <- sapply(1:p, function(idx) {
    if(idx %in% de_elements) {
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
  theta <- c(theta, c(rbind(beta.0, beta.1)))

  # size_factors.abundances <- sample(rep(GTEx_stats$total_read_profiles, ceiling(length(GTEx_stats$total_read_profiles)/(n*2))))[1:(n*2)]
  abundances <- sapply(1:p, function(j) {
    offset <- 2 + (j-1)*2 + 1
    mu.ij <- exp(theta[offset] + groups*theta[offset+1])
    rnbinom((n*2), mu = mu.ij, size = theta[1])
  })
  # samples x genes

  # transform these to proportions
  sampled_proportions <- apply(abundances, 1, function(x) x/sum(x))
  # genes x samples

  realized_total_counts <- rowSums(abundances)
  size_factors.observed_counts <- round(simcor(realized_total_counts, mean(realized_total_counts), sd(realized_total_counts),
                                               size_factor_correlation))
  # reflect any negative counts; this increases correlation but barely
  size_factors.observed_counts <- abs(size_factors.observed_counts)
  
  observed_counts <- sapply(1:(n*2), function(i) {
    if(is.na(size_factors.observed_counts[i])) {
      print(i)
    }
    rmultinom(1, size = size_factors.observed_counts[i], prob = sampled_proportions[,i])
  })
  # genes x samples

  return(list(abundances = abundances, observed_counts = t(observed_counts), groups = groups, theta = theta, de_genes = de_elements))
}

#' Simulate counts for control vs. treatment groups
#'
#' @param p number of genes
#' @param n number of samples per group
#' @param proportions set of proportions from which to sample true abundances
#' @param perturbed_features number of features to differentially express
#' @param fold_change fold change to induce in perturbed features
#' @param mean_size_factor average total abundance to simulate
#' @param size_factor_correlation correlation of observed abundances to original total abundances
#' @return named list of true abundances, observed abundances, and group assignments
#' @import entropy
#' @export
simulate_data <- function(p, n, proportions, perturbed_features, fold_change = NULL,
                          mean_size_factor = 5000, size_factor_correlation = 0) {
  groups <- c(rep(0, n), rep(1, n))
  theta <- c(0.2, 0.2) # fix the dispersions between groups for simplicity now
  beta.0 <- log(proportions)
  perturbed_elements <- sort(sample(1:p)[1:perturbed_features])
  if(!is.null(fold_change)) {
    scales <- rep(1, p)
    scales[perturbed_elements] <- fold_change
    beta.1 <- log(scales)
  } else {
    scales <- rep(0, p)
    random_signs <- rbinom(perturbed_features, size = 1, pro = rep(1/2, 2))
    random_signs[random_signs == 0] <- -1
    scales[perturbed_elements] <- random_signs*rnorm(perturbed_features, 1, 0.5)
    beta.1 <- scales
  }
  theta <- c(theta, c(rbind(beta.0, beta.1)))

  abundances <- matrix(NA, n*2, p)
  lower_bound <- mean_size_factor - mean_size_factor / 2
  upper_bound <- mean_size_factor + mean_size_factor / 2
  # scale size factor such that our beta.0 scale doesn't affect the total number of counts, i.e. if the mu.ij are
  #   small here (because we're sampling betas in large negative log space), increase the size factors to compensate
  #   such that we get total abundances at the same(ish) scale for all experiments
  # lower_bound <- lower_bound * (1/sum(exp(beta.0)))
  # upper_bound <- upper_bound * (1/sum(exp(beta.0)))
  size_factors.abundances <- round(runif(n*2, min = lower_bound, max = upper_bound))
  entropy_vector <- c()
  for(i in 1:(n*2)) {
    for(j in 1:p) {
      offset <- 2 + (j-1)*2 + 1
      mu.ij <- exp(theta[offset] + groups[i]*theta[offset+1])
      # if group dispersions are different...
      if(groups[i] == 0) {
        abundances[i,j] <- rnbinom(1, mu = size_factors.abundances[i]*mu.ij, size = 1 / theta[1])
      } else {
        abundances[i,j] <- rnbinom(1, mu = size_factors.abundances[i]*mu.ij, size = 1 / theta[2])
      }
    }
    entropy_vector <- c(entropy_vector, entropy(abundances[i,]))
  }

  mean_entropy <- mean(entropy_vector)
  max_entropy <- entropy(rep(1, p))
  relative_entropy <- mean_entropy / max_entropy

  # transform these to proportions
  sampled_proportions <- apply(abundances, 1, function(x) x/sum(x))

  # total observed counts
  realized_total_counts <- rowSums(abundances)
  size_factors.observed_counts <- round(simcor(realized_total_counts, mean(realized_total_counts), sd(realized_total_counts),
                                         size_factor_correlation))
  # reflect any negative counts; this increases correlation but barely
  size_factors.observed_counts <- abs(size_factors.observed_counts)
  # previously: wipe out all total abundance relationship
  # size_factors.observed_counts <- rowSums(abundances)[sample(1:(n*2))]

  observed_counts <- matrix(NA, n*2, p)
  for(i in 1:(n*2)) {
    observed_counts[i,] <- rmultinom(1, size = size_factors.observed_counts[i], prob = sampled_proportions[,i])
  }

  return(list(abundances = abundances, observed_counts = observed_counts, groups = groups,
              theta = theta, fold_change = fold_change, perturbed_features = perturbed_elements,
              size_factors.abundances = realized_total_counts, size_factors.observed_counts = size_factors.observed_counts,
              relative_entropy = relative_entropy))
}

#' Fit negative binomial model to null and full models and evaluate differential expression for a focal gene
#'
#' @param data simulated data set
#' @param feature_idx index of gene to test for differential expression
#' @param call_abundances if TRUE, call DE on abundances; if FALSE, on observed counts
#' @param rarefy if non-zero, this specifies the total abundance at which to resample the counts
#' @return p-value from NB GLM fit with MASS::glm.nb
#' @import MASS
#' @export
call_DE <- function(data, feature_idx, call_abundances = TRUE, rarefy = 0) {
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
    return(gene_data)
  }
}

#' Fit NBID model to null and full models and evaluate differential expression for a focal gene
#'
#' @param data simulated data set
#' @param feature_idx index of gene to test for differential expression
#' @param call_abundances if TRUE, call DE on abundances; if FALSE, on observed counts
#' @return p-value from likelihood ratio test of differential expression
#' @export
call_DE_original <- function(data, feature_idx, call_abundances = TRUE) {
  if(call_abundances) {
    base_params <- list(x = data$abundances[,feature_idx], groups = data$groups, size_factor = data$size_factors.abundances)
  } else {
    base_params <- list(x = data$observed_counts[,feature_idx], groups = data$groups, size_factor = data$size_factors.observed_counts)
  }
  
  # fit null model
  params_initial_H0 <- c(runif(2, min = 0.001, max = 1), runif(1, min = -1, max = 1))
  params_H0 <- base_params
  params_H0$null_model <- TRUE
  res_H0 <- optimize.NBID(params_initial_H0, params_H0)
  
  # fit full model
  params_initial_H1 <- c(params_initial_H0, runif(1, min = -1, max = 1))
  params_H1 <- base_params
  params_H1$null_model <- FALSE
  res_H1 <- optimize.NBID(params_initial_H1, params_H1)
  
  # apply LRT to give p-value
  return(calc_LRT.NBID(res_H0, res_H1, base_params))
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

#' Generate simulated data with randomly perturbed features, evaluate false negatives and positives in
#' differential expression calls and write these to an output file
#'
#' @param p number of genes
#' @param fold_change fold change to induce in perturbed features
#' @param proportion_perturbed_features proportion of features to differentially express
#' @param run_label human-readable string identifier for this simulation run
#' @param mean_size_factor average sequencing depth (total counts per sample)
#' @param size_factor_correlation correlation of observed abundances to original total abundances
#' @param bimodal if TRUE, simulates a composition with a high and low abundance cohorts
#' @param output_file output file to append to
#' @param alpha significant level below which to call a feature differentially expressed
#' @details Writes out error and simulation run statistics to a file.
#' @return NULL
#' @import LaplacesDemon
#' @import filelock
#' @export
run_evaluation_instance <- function(p, fold_change, proportion_perturbed_features, run_label, mean_size_factor = 5000,
                                    size_factor_correlation = 0, bimodal = FALSE, output_file = "results.txt", alpha = 0.05) {
  perturbed_features <- round(p*proportion_perturbed_features)
  all_taxa_observed <- FALSE
  while(!all_taxa_observed) {
    if(bimodal) {
      proportions <- rdirichlet(1, c(rep(25, round(p*0.2)), rep(1, round(p*0.8))))[1,]
    } else {
      proportions <- rdirichlet(1, rep(1, p))[1,]
    }
    data <- simulate_data(p = p, n = 50, proportions = proportions, perturbed_features = perturbed_features,
                          fold_change = fold_change, mean_size_factor = mean_size_factor, size_factor_correlation = size_factor_correlation)
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

  # calculate and compare differential expression calls
  calls.abundances <- c()
  calls.observed_counts <- c()
  for(i in 1:p) {
    pval.abundances <- call_DE(data, i, call_abundances = TRUE)
    pval.observed_counts <- call_DE(data, i, call_abundances = FALSE)
    if(!is.na(pval.abundances) & !is.na(pval.observed_counts)) {
      calls.abundances <- c(calls.abundances, pval.abundances <= alpha)
      calls.observed_counts <- c(calls.observed_counts, pval.observed_counts <= alpha)
    }
  }
  FN <- sum(calls.abundances == TRUE & calls.observed_counts == FALSE)
  FP <- sum(calls.abundances == FALSE & calls.observed_counts == TRUE)
  TN <- sum(calls.abundances == FALSE & calls.observed_counts == FALSE)
  TP <- sum(calls.abundances == TRUE & calls.observed_counts == TRUE)
  out_str <- paste0(run_label,"\t",
                    p,"\t",
                    length(data$perturbed_features),"\t",
                    round(data$relative_entropy, 5),"\t",
                    fold_change,"\t",
                    mean_size_factor,"\t",
                    size_factor_correlation,"\t",
                    bimodal,"\t",
                    alpha,"\t",
                    TP,"\t",
                    FP,"\t",
                    TN,"\t",
                    FN)
  record_result(out_str, file.path("output","results.tsv"))
}

#' Generate simulated data with RNA-seq like features and > 3 fold differential expression, evaluate false negatives and positives in
#' differential expression calls and write these to an output file
#'
#' @param p number of genes
#' @param n number of samples per group
#' @param proportion_de proportion of differentially expressed genes
#' @param run_label human-readable string identifier for this simulation run
#' @param size_factor_correlation correlation of observed abundances to original total abundances
#' @param output_file output file to append to
#' @param alpha significant level below which to call a feature differentially expressed
#' @details Writes out error and simulation run statistics to a file.
#' @return NULL
#' @export
run_RNAseq_evaluation_instance <- function(p, n, proportion_de, run_label, size_factor_correlation = 0, output_file = "results.txt", alpha = 0.05) {
  all_taxa_observed <- FALSE
  while(!all_taxa_observed) {
    cat("Simulating data...\n")
    data <- simulate_RNAseq(p = p, n = n, proportion_de = proportion_de, size_factor_correlation = size_factor_correlation)
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

  # calculate and compare differential expression calls
  cat("Evaluating differential expression...\n")
  calls.abundances <- c()
  calls.observed_counts <- c()
  for(i in 1:p) {
    pval.abundances <- call_DE(data, i, call_abundances = TRUE)
    pval.observed_counts <- call_DE(data, i, call_abundances = FALSE)
    if(!is.na(pval.abundances) & !is.na(pval.observed_counts)) {
      calls.abundances <- c(calls.abundances, pval.abundances <= alpha)
      calls.observed_counts <- c(calls.observed_counts, pval.observed_counts <= alpha)
    }
  }
  FN <- sum(calls.abundances == TRUE & calls.observed_counts == FALSE)
  FP <- sum(calls.abundances == FALSE & calls.observed_counts == TRUE)
  TN <- sum(calls.abundances == FALSE & calls.observed_counts == FALSE)
  TP <- sum(calls.abundances == TRUE & calls.observed_counts == TRUE)
  out_str <- paste0(run_label,"\t",
                    p,"\t",
                    proportion_de,"\t",
                    size_factor_correlation,"\t",
                    TP,"\t",
                    FP,"\t",
                    TN,"\t",
                    FN)
  record_result(out_str, file.path("output","results_RNAseq.tsv"))
}

#' Generate simulated data with RNA-seq like features and > 3 fold differential expression, evaluate false negatives and positives in
#' differential expression calls and write these to an output file
#'
#' @param p number of genes
#' @param n number of samples per group
#' @param run_label human-readable string identifier for this simulation run
#' @param de_sweep vector of proportions of differentially expressed genes to simulate
#' @param corr_sweep vector of correlations of size factors (true vs. observed abundances) to simulate
#' @param output_file output file to append to
#' @param alpha significant level below which to call a feature differentially expressed
#' @details Writes out error and simulation run statistics to a file.
#' @return NULL
#' @export
sweep_RNAseq <- function(p, n, run_label, de_sweep = seq(from = 0.1, to = 0.9, by = 0.1),
                         corr_sweep = seq(from = 0.1, to = 0.9, by = 0.1), output_file = "results.txt", alpha = 0.05) {
  for(de_prop in de_sweep) {
    for(sf_corr in corr_sweep) {
      cat("Evaluating size factor correlation =",round(sf_corr, 2),"and DE proportion =",round(de_prop, 2),"\n")
      run_RNAseq_evaluation_instance(p, n, proportion_de = de_prop, run_label = "RNAseq_like_sweep", size_factor_correlation = sf_corr,
                                     output_file = "results.txt", alpha = 0.05)
    }
  }
}









