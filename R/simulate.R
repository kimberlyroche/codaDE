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
#' @param library_size_correlation correlation of observed abundances to original absolute abundances
#' @param spike_in flag indicating whether or not to simulate a spike-in (control) with very low dispersion
#' @return named list of true abundances and observed abundances
#' @export
resample_counts <- function(theta, groups, library_size_correlation, spike_in = FALSE) {
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
                                               library_size_correlation))
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
#' @param n number of samples per group
#' @param p number of features (i.e. genes or bacterial taxa)
#' @param k number of distinct cell types
#' @param ref_data specifies the reference data set to use ("Morton", "Athanasiadou_ciona", "simulated"; case-sensitive)
#' @param asymmetry proportion of baselines to draw from condition 1
#' @param sequencing_depth average total nUMI to simulate
#' @param proportion_da proportion of differentially abundant genes
#' @param library_size_correlation correlation of observed abundances to original total abundances (see details)
#' @param spike_in simulate a control spike in with low dispersion
#' @param possible_fold_changes if specified, a distribution of possible fold changes available for differential abundance simulation
#' @return named list of true abundances, observed abundances, and group assignments
#' @details Library size correlation can be specified as 0, 1, or 2, indicating no correlation, modest correlation, or very strong
#' correlation with true abundances.
#' 
#' Samples the mean gene abundances from Halpern et al. (2017) as a references. This gives a data set
#' with on average 75% zeros which may still not be representative in terms of sparsity. Need to follow up. (7/22/2020)
#' @import LaplacesDemon
#' @export
simulate_sequence_counts <- function(n = 500, p = 1000, k = 1, ref_data = "Athanasidou_ciona", asymmetry = 0.5,
                                     sequencing_depth = 1e5, proportion_da = 0.1, library_size_correlation = 0,
                                     spike_in = FALSE, possible_fold_changes = NULL) {

  # We'll sample out of GTEx protein-coding gene estimates
  # GTEx_estimates_baseline <- readRDS("data/empirical_mu_size.rds")
  #
  # if(microbial) {
  #   estimate_pairs <- readRDS("data/empirical_mu_size_pairs_16S.rds")
  #   estimate_pairs <- cbind(estimate_pairs[,3:4], estimate_pairs[,1:2])
  #   k <- 1
  #   p <- 1000 # SV count
  # } else {
  #   estimate_pairs <- readRDS("data/empirical_mu_size_pairs.rds")
  #   p <- 20000 # gene count
  # }
  # # Change zero-mean to 1, still rare
  # estimate_pairs[estimate_pairs[,1] == 0,1] <- min(estimate_pairs[estimate_pairs[,1] > 0,1])
  # estimate_pairs[estimate_pairs[,3] == 0,3] <- min(estimate_pairs[estimate_pairs[,3] > 0,3])
  # # Fix unestimable (very low) dispersions
  # estimate_pairs[is.na(estimate_pairs[,2]),2] <- 100
  # estimate_pairs[is.na(estimate_pairs[,4]),4] <- 100
  # 
  # n_da <- round(p * proportion_da)
  # n_not_da <- p - n_da
  # if(!is.null(possible_fold_changes)) {
  #   da_assignment <- as.logical(rbinom(p, size = 1, prob = proportion_da))
  #   # Randomly draw mean-dispersions from across tissues in the empirical sample for non-differential genes
  #   all_params <- unname(estimate_pairs[sample(1:nrow(estimate_pairs), size = p, replace = TRUE),])
  #   all_params <- cbind(all_params, all_params) # before-after are identical but...
  #   da_factors <- sample(c(1/possible_fold_changes, possible_fold_changes), size = p, replace = TRUE)
  #   for(pp in 1:p) {
  #     if(da_assignment[pp]) {
  #       all_params[pp,3] <- all_params[pp,3]*da_factors[pp]
  #       # Dispersion unchanged
  #     }
  #   }
  # } else {
  #   da_assignment <- c(rep(FALSE, n_not_da), rep(TRUE, n_da))
  #   # Randomly draw mean-dispersions from across tissues in the empirical sample for non-differential genes
  #   non_da_params <- unname(estimate_pairs[sample(1:nrow(estimate_pairs), size = n_not_da, replace = TRUE),1:2])
  #   non_da_params <- cbind(non_da_params, non_da_params) # before-after are identical
  #   # Randomly draw in a similar way for the differential set
  #   # Require a fold change in means that's < 1/2 or > 2
  #   ratio_means <- estimate_pairs[,3] / estimate_pairs[,1]
  #   thresholded_genes <- which(ratio_means <= 0.5 | ratio_means >= 2)
  #   da_params <- unname(estimate_pairs[sample(thresholded_genes, size = n_da, replace = TRUE),])
  #   all_params <- rbind(non_da_params, da_params)
  # }
  #
  # Clean up NAs
  # all_params[is.na(all_params[,2]),2] <- 1
  # all_params[is.na(all_params[,4]),4] <- 1
  
  ref_file <- file.path("data", paste0("DE_reference_",ref_data,".rds"))
  cat("Using reference file:",ref_file,"\n")
  data_obj <- readRDS(ref_file)
  # counts <- data_obj$counts
  # baseline_distribution <- data_obj$baseline_distribution
  # differential_distribution <- data_obj$differential_distribution
  
  n_da <- round(p * proportion_da)
  n_not_da <- p - n_da

  if(!is.null(possible_fold_changes)) {
    # I think this is broken -- worth a look
    # sampled_features <- sample(1:length(data_obj$cond1), size = p, replace = TRUE)
    # all_params <- matrix(NA, p, 2)
    # for(i in 1:length(sampled_features)) {
    #   all_params[i,] <- c(data_obj$cond1[sampled_features[i]], 1)
    # }
    # all_params <- cbind(all_params, all_params) # before-after are identical but...
    # 
    # da_assignment <- as.logical(rbinom(p, size = 1, prob = proportion_da))
    # # Randomly draw mean-dispersions from across tissues in the empirical sample for non-differential genes
    # da_factors <- sample(c(1/possible_fold_changes, possible_fold_changes), size = p, replace = TRUE)
    # for(pp in 1:p) {
    #   if(da_assignment[pp]) {
    #     all_params[pp,3] <- all_params[pp,3]*da_factors[pp]
    #     # Dispersion unchanged
    #   }
    # }
  } else {
    # differential_baselines <- data_obj$differential_distribution[,1]
    # breaks <- quantile(differential_baselines, probs = seq(from = 0, to = 1, length.out = 5))
    # # Fix for R's cut function being kind of shitty
    # breaks[[1]] <- -Inf
    # breaks[[length(breaks)]] <- Inf
    # 
    # da_sample <- sample(1:nrow(all_params), size = n_da)
    # da_assignment <- logical(p)
    # da_assignment[da_sample] <- TRUE
    # 
    # for(baseline_idx in da_sample) {
    #   baseline_mu <- all_params[baseline_idx,1]
    #   baseline_size <- all_params[baseline_idx,2]
    #   
    #   binned_expression <- cut(c(baseline_mu, differential_baselines), breaks = breaks)
    #   baseline_mu_binned <- binned_expression[1]
    #   binned_expression <- binned_expression[2:length(binned_expression)]
    #   similar_mus <- which(baseline_mu_binned == binned_expression)
    #   
    #   differential_idx <- sample(similar_mus, size = 1)
    #   differential_mu <- differential_distribution[differential_idx,3]
    #   differential_size <- differential_distribution[differential_idx,4]
    #   
    #   all_params[baseline_idx,3:4] <- c(differential_mu, differential_size)
    #   
    #   # Visually validate a "difference"
    #   # expr1 <- rnbinom(100, mu = baseline_mu, size = baseline_size)
    #   # expr2 <- rnbinom(100, mu = differential_mu, size = differential_size)
    #   # plot(c(expr1, expr2))
    # }
    
    # Assume symmetric sampling
    # if(p %% 2 != 0) { p <- p + 1 }
    
    p1 <- round(asymmetry*p)
    p2 <- p - p1
    
    # baseline_features <- sample(1:length(data_obj$cond1), size = p, replace = TRUE)
    # baseline_features1 <- baseline_features[1:(p/2)]
    # baseline_features2 <- baseline_features[(p/2 + 1):p]
    baseline_features1 <- sample(1:length(data_obj$cond1), size = p1, replace = TRUE)
    baseline_features2 <- sample(1:length(data_obj$cond2), size = p2, replace = TRUE)
    baseline_counts <- c(data_obj$cond1[baseline_features1], data_obj$cond2[baseline_features2])
    treatment_counts <- c(data_obj$cond2[baseline_features1], data_obj$cond1[baseline_features2])
    
    da_sample <- sample(1:p, size = n_da)
    da_assignment <- logical(p)
    da_assignment[da_sample] <- TRUE
    
    all_params <- matrix(1, p, 4)
    all_params[,c(1,3)] <- baseline_counts
    all_params[,c(2,4)] <- 1000 # fixed dispersion for now
    all_params[da_assignment,3] <- treatment_counts[da_assignment]
    
    # Visually validate a "difference"
    # idx <- sample(da_sample, size = 1)
    # expr1 <- rnbinom(100, mu = all_params[idx,1], size = all_params[idx,2])
    # expr2 <- rnbinom(100, mu = all_params[idx,3], size = all_params[idx,4])
    # plot(c(expr1, expr2))
  }
  
  # Assign control vs. treatment groups
  groups <- c(rep(0, n), rep(1, n))
  # groups <- as.factor(da_assignment)
  # levels(groups) <- c(0, 1)
  
  # Sample true counts (format: samples [rows] x genes [columns])
  abundances <- sapply(1:p, function(j) {
    counts <- rnbinom(n*2,
                      mu = c(rep(all_params[j,1], n), rep(all_params[j,3], n)),
                      size = c(rep(all_params[j,2], n), rep(all_params[j,4], n)))
    counts
  })

  realized_total_counts <- rowSums(abundances)
  scale_down_factor <- sequencing_depth / mean(realized_total_counts)
  scaled_real_counts <- round(scale_down_factor*realized_total_counts)
  
  # Rescale the observed counts in advance of proportional resampling
  if(library_size_correlation == 2) {
    # Very strong correlation
    library_sizes.observed_counts <- rnbinom(n*2, mu = scaled_real_counts, size = 100)
  } else if(library_size_correlation == 1) {
    # Modest correlation
    library_sizes.observed_counts <- rnbinom(n*2, mu = scaled_real_counts, size = 10)
  } else {
    # No correlation
    library_sizes.observed_counts <- rnbinom(n*2, mu = sequencing_depth, size = 100)
  }
  
  # transform to proportions
  sampled_proportions <- apply(abundances, 1, function(x) x/sum(x))
  observed_counts <- sapply(1:(n*2), function(i) {
    resample_probs <- sampled_proportions[,i]
    # censor very low concentration stuff to simulate drop-outs
    # I'm trying to simulate some minimum detectable concentration
    # this threshold is pretty arbitrary; we need to tune this to resemble real data
    resample_probs[resample_probs < 1e-8] <- 0
    rmultinom(1, size = library_sizes.observed_counts[i], prob = resample_probs)
  })
  observed_counts <- t(observed_counts)
  
  # Quick visualization
  # par(mfrow = c(1, 2))
  # idx <- sample(which(da_assignment == FALSE))[1]
  # plot(1:n, abundances[1:n,idx], xlim = c(0, 2*n), ylim = c(0, max(abundances[,idx])))
  # lines((n+1):(2*n), abundances[(n+1):(2*n),idx], type = "p", col = "red")
  # plot(1:n, observed_counts[1:n,idx], xlim = c(0, 2*n), ylim = c(0, max(observed_counts[,idx])))
  # lines((n+1):(2*n), observed_counts[(n+1):(2*n),idx], type = "p", col = "red")

  return(list(abundances = abundances,
              observed_counts = observed_counts,
              groups = groups,
              da_assignment = which(da_assignment == TRUE)))
}
