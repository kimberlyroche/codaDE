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

# Utility function to resample counts
resample_counts <- function(proportions, library_sizes) {
  observed_counts <- sapply(1:ncol(proportions), function(i) {
    resample_probs <- proportions[,i]
    # censor very low concentration stuff to simulate drop-outs
    # I'm trying to simulate some minimum detectable concentration
    # this threshold is pretty arbitrary; we need to tune this to resemble real data
    resample_probs[resample_probs < 1e-8] <- 0
    rmultinom(1, size = library_sizes[i], prob = resample_probs)
  })
  t(observed_counts)
}

# Utility function to spike-in a one-count to prevent any fully unobserved taxa
# Counts should be samples x features
spike_in_ones <- function(counts) {
  all_zeros <- which(colSums(counts) == 0)
  n <- nrow(counts)/2
  for(idx in all_zeros) {
    one_idx <- sample(1:n, size = 1)
    counts[one_idx,idx] <- 1
    one_idx <- sample((n+1):(n*2), size = 1)
    counts[one_idx,idx] <- 1
  }
  return(counts)
}

#' Simulate UMI count single-cell RNA-seq data
#'
#' @param n number of samples per group
#' @param p number of features (i.e. genes or bacterial taxa)
#' @param ref_data_name optional: specifies the reference data set to use
#' ("Morton", "Athanasiadou_ciona", "simulated"; case-sensitive)
#' @param data_obj reference data set object; if non-null, this takes precedence
#' over ref_data_name
#' @param asymmetry proportion of baselines to draw from condition 1
#' @param proportion_da proportion of differentially abundant genes
#' @param spike_in simulate a control spike in with low dispersion
#' @return named list of true abundances, observed abundances, and group assignments
#' @details Library size correlation can be specified as 0, 1, or 2, indicating no correlation, modest correlation, or very strong
#' correlation with true abundances.
#' 
#' Samples the mean gene abundances from Halpern et al. (2017) as a references. This gives a data set
#' with on average 75% zeros which may still not be representative in terms of sparsity. Need to follow up. (7/22/2020)
#' @import LaplacesDemon
#' @export
simulate_sequence_counts <- function(n = 500,
                                     p = 1000,
                                     ref_data_name = "Athanasidou_ciona",
                                     data_obj = NULL,
                                     asymmetry = 0.5,
                                     proportion_da = 0.1,
                                     spike_in = FALSE) {
  if(is.null(ref_data_name) & is.null(data_obj)) {
    stop("One of ref_data_name and ref_data must be specified!")
  }
  if(is.null(data_obj)) {
    ref_file <- file.path("data", paste0("DE_reference_",ref_data_name,".rds"))
    cat("Using reference file:",ref_file,"\n")
    data_obj <- readRDS(ref_file)
  }
  # counts <- data_obj$counts
  # baseline_distribution <- data_obj$baseline_distribution
  # differential_distribution <- data_obj$differential_distribution
  
  n_da <- round(p * proportion_da)
  n_not_da <- p - n_da

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

  if(spike_in) {
    spike_in_mean <- mean(all_params[,1]) # mean baseline abundance
    # No noise
    # abundances <- cbind(abundances, rep(spike_in_mean, n*2))
    # Noisy
    abundances <- cbind(abundances, rnbinom(n*2, mu = spike_in_mean, size = 100))
  }

  # If there are any fully-zero taxa, spike-in a random 1-count in each condition
  # This is so these taxa don't get summarily filtered out by some of the DA
  # calling methods. It basically helps with matching p-value output to indices.
  abundances <- spike_in_ones(abundances)
  
  realized_total_counts <- rowSums(abundances)
  # scale_down_factor <- sequencing_depth / mean(realized_total_counts)
  # scaled_real_counts <- round(scale_down_factor*realized_total_counts)
  
  # On average, no correlation in true / observed library sizes
  # library_sizes.observed_counts1 <- rnbinom(n*2, mu = realized_total_counts, size = 0.1)
  # Modest correlation true / observed library sizes
  # library_sizes.observed_counts2 <- rnbinom(n*2, mu = realized_total_counts, size = 50)
  # # Very strong correlation in true / observed library sizes
  # library_sizes.observed_counts3 <- rnbinom(n*2, mu = realized_total_counts, size = 100)

  library_sizes.observed_counts1 <- sample(realized_total_counts)
  # library_sizes.observed_counts2 <- realized_total_counts
  # shuffle_idx <- sample(1:(n*2), size = n)
  # library_sizes.observed_counts2[shuffle_idx] <- library_sizes.observed_counts2[sample(shuffle_idx)]
  # library_sizes.observed_counts2 <- abs(rnorm(n*2, mean = realized_total_counts, sd = sd(realized_total_counts)*2))
  
  # Let's try a little less variation
  library_sizes.observed_counts2 <- rnorm(n*2, mean = realized_total_counts, sd = sd(realized_total_counts))
  library_sizes.observed_counts2 <- sapply(library_sizes.observed_counts2, function(x) max(x, 1000))
  
  # ggplot(data = data.frame(x = rep(1:(n*2), 2),
  #                          totals = c(realized_total_counts, library_sizes.observed_counts2),
  #                          type = c(rep("A", n*2), rep("B", n*2))),
  #        aes(x = x, y = totals, color = type)) +
  #   geom_point(size = 3, alpha = 0.5)
  # 
  # cor(realized_total_counts, library_sizes.observed_counts2)
  
  # transform to proportions
  sampled_proportions <- apply(abundances, 1, function(x) x/sum(x))
  observed_counts1 <- resample_counts(sampled_proportions, library_sizes.observed_counts1)
  observed_counts2 <- resample_counts(sampled_proportions, library_sizes.observed_counts2)
  # observed_counts3 <- resample_counts(sampled_proportions, library_sizes.observed_counts3)
  
  # Repeat spike-in of ones
  observed_counts1 <- spike_in_ones(observed_counts1)
  observed_counts2 <- spike_in_ones(observed_counts2)
  
  # Quick visualization
  # par(mfrow = c(1, 2))
  # idx <- sample(which(da_assignment == FALSE))[1]
  # plot(1:n, abundances[1:n,idx], xlim = c(0, 2*n), ylim = c(0, max(abundances[,idx])))
  # lines((n+1):(2*n), abundances[(n+1):(2*n),idx], type = "p", col = "red")
  # plot(1:n, observed_counts[1:n,idx], xlim = c(0, 2*n), ylim = c(0, max(observed_counts[,idx])))
  # lines((n+1):(2*n), observed_counts[(n+1):(2*n),idx], type = "p", col = "red")

  return(list(spike_in = spike_in,
              baseline_counts = baseline_counts,
              treatment_counts = treatment_counts,
              abundances = abundances,
              observed_counts1 = observed_counts1,
              observed_counts2 = observed_counts2,
              # observed_counts3 = observed_counts3,
              groups = groups,
              da_assignment = which(da_assignment == TRUE)))
}
