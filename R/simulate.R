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
#' @param asymmetry proportion of baseline abundances to draw from condition 1
#' @param proportion_da proportion of differentially abundant genes
#' @param spike_in simulate a control feature with low dispersion across 
#' conditions; this is indexed as the last feature
#' @return named list of true abundances, observed abundances, and group 
#' assignments
#' @import LaplacesDemon
#' @export
simulate_sequence_counts <- function(n = 20,
                                     p = 100,
                                     ref_data_name = "simulated",
                                     data_obj = NULL,
                                     asymmetry = 1,
                                     proportion_da = 0.75,
                                     spike_in = FALSE) {
  if(is.null(ref_data_name) & is.null(data_obj)) {
    stop("One of ref_data_name and ref_data must be specified!")
  }
  if(is.null(data_obj)) {
    ref_file <- file.path("data", paste0("DE_reference_",ref_data_name,".rds"))
    data_obj <- readRDS(ref_file)
  }

  n_da <- round(p * proportion_da)
  n_not_da <- p - n_da

  p1 <- round(asymmetry*p)
  p2 <- p - p1
  
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
  
  # Sample true counts (format: samples [rows] x genes [columns])
  abundances <- sapply(1:p, function(j) {
    counts <- rnbinom(n*2,
                      mu = c(rep(all_params[j,1], n), rep(all_params[j,3], n)),
                      size = c(rep(all_params[j,2], n), rep(all_params[j,4], n)))
    counts
  })

  if(spike_in) {
    spike_in_mean <- mean(all_params[,1]) # mean baseline abundance
    # Noisy
    abundances <- cbind(abundances, rnbinom(n*2, mu = spike_in_mean, size = 100))
  }

  # If there are any fully-zero taxa, spike-in a random 1-count in each 
  # condition. This is so these taxa don't get summarily filtered out by some of 
  # the DA calling methods. It basically helps with matching p-value output to 
  # indices.
  abundances <- spike_in_ones(abundances)
  
  realized_total_counts <- rowSums(abundances)
  
  # These are relative counts
  library_sizes.observed_counts1 <- sample(realized_total_counts)
  
  # These are counts with total abundances correlated with true abundances
  library_sizes.observed_counts2 <- rnorm(n*2,
                                          mean = realized_total_counts,
                                          sd = sd(realized_total_counts))
  library_sizes.observed_counts2 <- sapply(library_sizes.observed_counts2,
                                           function(x) max(x, 1000))
  library_sizes.observed_counts2 <- round(library_sizes.observed_counts2)
  
  # Transform to proportions (samples x features)
  sampled_proportions <- t(apply(abundances, 1, function(x) x/sum(x)))
  observed_counts1 <- resample_counts(sampled_proportions,
                                      library_sizes.observed_counts1)
  observed_counts2 <- resample_counts(sampled_proportions,
                                      library_sizes.observed_counts2)

  # Repeat spike-in of ones
  observed_counts1 <- spike_in_ones(observed_counts1)
  observed_counts2 <- spike_in_ones(observed_counts2)
  
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
