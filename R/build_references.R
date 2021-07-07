#' Generate simulated differential expression for two conditions
#'
#' @param p number of features (genes, taxa) to simulate
#' @param log_mean log mean for condition 1
#' @param log_var log variance for condition 1
#' @param log_noise_var log variance for noise component added to condition 1
#' to give condition 2
#' @param base_correlation this is the base correlation matrix for simulated
#' features in log space
#' @param concentration concentration parameter for sampling the correlation 
#' across feature; if base_correlation is NULL and concentration = LARGE, features
#' are effectively independent
#' @param save_name optional tag to append to saved data file
#' @return NULL
#' @import matrixsampling
#' @import MASS
#' @export
build_simulated_reference <- function(p = 1000, log_mean = 0, log_var = 4,
                                      log_noise_var = 4, base_correlation = NULL,
                                      concentration = 1e6, save_name = NULL) {
  if(is.null(base_correlation)) {
    base_correlation <- diag(p)
  }
  if(concentration < p + 2) {
    stop("Invalid concentration specified!")
  }
  # log_counts1 <- rnorm(p, log_mean, log_var)
  # log_counts2 <- log_counts1 + rnorm(p, 0, log_noise_var)
  K <- cov2cor(rinvwishart(1, concentration, base_correlation*concentration)[,,1])

  # Add randomness to log perturbation size
  # log_noise_var2 <- runif(1, min = 0.1, max = log_noise_var)
  # log_counts1 <- mvrnorm(1, rep(log_mean, p), diag(p)*log_var)
  log_counts1 <- rnorm(p, rep(log_mean, p), sqrt(log_var))
  log_perturbation <- mvrnorm(1, rep(0, p), K*log_noise_var)
  log_counts2 <- log_counts1 + log_perturbation

  counts1 <- exp(log_counts1)
  counts2 <- exp(log_counts2)
  
  sim_obj <- list(log_cond1 = log_counts1, log_cond2 = log_counts2,
                  cond1 = counts1, cond2 = counts2,
                  log_perturbation = log_perturbation, correlation_matrix = K)
  
  if(is.null(save_name)) {
    return(sim_obj)
  } else {
    # save_file <- ifelse(is.null(save_name),
    #                     "DE_reference_simulated.rds",
    #                     paste0("DE_reference_",save_name,".rds"))
    save_file <- paste0("DE_reference_",save_name,".rds")
    saveRDS(sim_obj,
            file = file.path("data", save_file))
  }
}

#' Generate differential expression reference for Barlow et al. (2020) 16S data
#'
#' @return NULL
#' @export
build_Barlow_reference <- function() {
  parsed_obj <- parse_Barlow()
  counts <- parsed_obj$counts
  groups <- parsed_obj$groups
  
  # Build DE model
  counts1 <- counts[,groups == "control" & days == 10]
  counts2 <- counts[,groups == "keto" & days == 10]
  
  # Arbitrarily "match" a subject on day 10 in condition 1 with a subject on
  # day 10 in condition 2
  
  counts1 <- c(counts1)
  counts2 <- c(counts2)
  
  saveRDS(list(cond1 = counts1, cond2 = counts2),
          file = file.path("data", "DE_reference_Barlow.rds"))
  return(NULL)
}

#' Generate differential expression reference for Morton et al. (2019) 16S data
#'
#' @return NULL
#' @export
build_Morton_reference <- function() {
  parsed_obj <- parse_Morton() # check this -- this fn. has been altered!
  counts <- parsed_obj$counts
  groups <- parsed_obj$groups

  # Build DE model
  counts1 <- c(counts[,groups == "before"])
  counts2 <- c(counts[,groups == "after"])

  saveRDS(list(cond1 = counts1, cond2 = counts2),
          file = file.path("data", "DE_reference_Morton.rds"))
  return(NULL)
}

#' Generate differential expression reference for Athanasiadou et al. (2021) 
#' bulk RNA-seq data
#'
#' @return named list of observed counts in conditions 1 and 2
#' @export
build_Athanasiadou_reference <- function() {
  # TBD - allow choice of baseline v. treatment
  parsed_obj <- parse_Athanasiadou(which_data = "yeast")
  counts <- parsed_obj$counts
  groups <- parsed_obj$groups

  # Build DE model (arbitrarily matching replicates)
  counts1 <- c(counts[,groups == "C12"])
  counts2 <- c(counts[,groups == "C30"])
  
  saveRDS(list(cond1 = counts1, cond2 = counts2),
          file = file.path("data", "DE_reference_Athanasiadou_yeast.rds"))

  parsed_obj <- parse_Athanasiadou(which_data = "ciona")
  counts <- parsed_obj$counts
  groups <- parsed_obj$groups
  
  # Build DE model
  counts1 <- c(counts[,groups == "lacz"])
  counts2 <- c(counts[,groups == "dnfgfr"])
  
  saveRDS(list(cond1 = counts1, cond2 = counts2),
          file = file.path("data", "DE_reference_Athanasiadou_ciona.rds"))
}
