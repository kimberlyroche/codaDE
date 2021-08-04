#' Utility function to multinomially resample counts
#'
#' @param proportions relative abundance table (samples x features)
#' @return resampled count table
#' @export
resample_counts <- function(proportions, library_sizes) {
  observed_counts <- sapply(1:nrow(proportions), function(i) {
    resample_probs <- proportions[i,]
    # Censor very low concentration stuff to simulate drop-outs. I'm trying to 
    # simulate some minimum detectable concentration but this threshold is 
    # pretty arbitrary. We need to tune this to resemble real data.
    resample_probs[resample_probs < 1e-8] <- 0
    rmultinom(1, size = library_sizes[i], prob = resample_probs)
  })
  t(observed_counts)
}

#' Utility function to spike-in a one-count to prevent any fully unobserved taxa
#'
#' @param counts count table (samples x features)
#' @return counts with a minimum single observation of each feature
#' @export
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
#' @param ref_data_name optional name of reference data set to use from among
#' "Morton", "Athanasiadou_ciona", "simulated"
#' @param data_obj optional reference data set object; this takes precedence 
#' over ref_data_name
#' @param replicate_noise optional deviation parameter associated with variation
#' in replicate totals (a max of 1 is approx. an upper limit)
#' @param proportion_da optional proportion of differentially abundant genes
#' @return named list of true abundances, observed abundances, and group 
#' assignments
#' @import LaplacesDemon
#' @export
simulate_sequence_counts <- function(n = 20,
                                     p = 100,
                                     ref_data_name = "simulated",
                                     data_obj = NULL,
                                     replicate_noise = 0,
                                     proportion_da = 0.75) {
  if(is.null(ref_data_name) & is.null(data_obj)) {
    stop("One of ref_data_name and ref_data must be specified!")
  }
  if(is.null(data_obj)) {
    ref_file <- file.path("data", paste0("DE_reference_",ref_data_name,".rds"))
    data_obj <- readRDS(ref_file)
  }

  n_da <- round(p * proportion_da)
  n_not_da <- p - n_da

  # Randomly select which condition to use as the baseline
  if(rbinom(1, size = 1, prob = 0.5) == 1) {
    # Base condition is #1
    p1 <- p
    p2 <- 0
  } else {
    # Base condition is #2
    p1 <- 0
    p2 <- p
  }
  
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
  
  # Assign control vs. treatment groups
  groups <- c(rep(0, n), rep(1, n))
  
  mu_scale <- sapply(rnorm(n*2, 1, replicate_noise), function(x) max(0.1, x))
  # Sample true counts (format: samples [rows] x genes [columns])
  abundances <- sapply(1:p, function(j) {
    counts <- rnbinom(n*2,
                      mu = c(rep(all_params[j,1], n), rep(all_params[j,3], n))*mu_scale,
                      size = c(rep(all_params[j,2], n), rep(all_params[j,4], n)))
    counts
  })

  # Spike-in; not used
  # spike_in_mean <- mean(all_params[,1]) # mean baseline abundance
  # abundances <- cbind(abundances, rnbinom(n*2, mu = spike_in_mean, size = 100))

  # If there are any fully-zero taxa, spike-in a random 1-count in each 
  # condition. This is so these taxa don't get summarily filtered out by some of 
  # the DA calling methods. It basically helps with matching p-value output to 
  # indices.
  abundances <- spike_in_ones(abundances)
  
  realized_total_counts <- rowSums(abundances)
  
  # Arbitrarily choose an average sequencing depth
  # A Poisson has awfully little variation relative to "real" data
  # library_sizes.observed_counts1 <- rpois(n*2, runif(1, min = 5000, 2e06))
  library_sizes.observed_counts1 <- rnbinom(n*2, mu = runif(1, min = 5000, 2e06), size = 100)
  
  # Choose an average sequencing depth centered on the original
  # library_sizes.observed_counts2 <- c(rnbinom(n, mu = mean(realized_total_counts[1:n]), size = 0.5),
  #                                     rnbinom(n, mu = mean(realized_total_counts[(n+1):(n*2)]), size = 0.5))
  # library_sizes.observed_counts2 <- sapply(library_sizes.observed_counts2, function(x) {
  #   if(x < 5000) {
  #     5000
  #   } else if(x > 2e06) {
  #     2e06
  #   } else {
  #     x
  #   }
  # })
  
  # Calculate the proportional differential
  FC <- (mean(realized_total_counts[(n+1):(n*2)]) - mean(realized_total_counts[1:n])) / mean(realized_total_counts)
  # percent_FC_retained <- runif(1, min = 0.2, max = 0.8)
  percent_FC_retained <- 0.9
  print(percent_FC_retained)
  sampled_FC <- FC*percent_FC_retained
  
  new_baseline_depth <- runif(1, min = 5000, 1e06)
  new_differential <- new_baseline_depth * sampled_FC
  
  if(new_differential > 0) {
    # B larger than A
    depth_A <- new_baseline_depth
    depth_B <- depth_A + new_differential
  } else {
    # B smaller than A
    depth_B <- new_baseline_depth
    depth_A <- depth_B + abs(new_differential)
  }
  
  library_sizes.observed_counts2 <- c(rnbinom(n, mu = depth_A, size = 10),
                                      rnbinom(n, mu = depth_B, size = 10))
  library_sizes.observed_counts2 <- sapply(library_sizes.observed_counts2, function(x) {
    if(x < 5000) {
      5000
    } else if(x > 2e06) {
      2e06
    } else {
      x
    }
  })
  
  # Transform to proportions (samples x features)
  sampled_proportions <- t(apply(abundances, 1, function(x) x/sum(x)))
  observed_counts1 <- resample_counts(sampled_proportions,
                                      library_sizes.observed_counts1)
  observed_counts2 <- resample_counts(sampled_proportions,
                                      library_sizes.observed_counts2)

  # Repeat spike-in of ones
  observed_counts1 <- spike_in_ones(observed_counts1)
  observed_counts2 <- spike_in_ones(observed_counts2)
  
  return(list(baseline_counts = baseline_counts,
              treatment_counts = treatment_counts,
              abundances = abundances,
              observed_counts1 = observed_counts1,
              observed_counts2 = observed_counts2,
              # observed_counts3 = observed_counts3,
              groups = groups,
              da_assignment = which(da_assignment == TRUE)))
}
