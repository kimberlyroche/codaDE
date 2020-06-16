#' Calculate negative log likelihood for NBID model
#'
#' @param theta vector of parameters (see details)
#' @param other_params named list of counts, group assignments, size factor (per-sample total counts), and H0 flag
#' @details Parameter vector contains (in this order) control group dispersion, treatment group dispersion,
#' focal gene intercept, and (if not evaluating null model) focal gene treatment group effect.
#' @return negative log likelihood for focal gene
#' @export
#' @examples
#' theta <- c(runif(2, min = 0.001, max = 1), runif(2, min = -1, max = 1))
#' other_params <- list(x = observed_counts, groups = sample_group_assignments, null_model = FALSE)
#' nll <- negativeLL.NBID(theta, other_params)
negativeLL.NBID <- function(theta, other_params) {
  x <- other_params$x
  groups <- other_params$group
  size_factor <- other_params$size_factor # this means mu is learned as a proportion
  #size_factor <- rep(mean(x), length(x)) # this means mu is learned as a scale
  null_model <- other_params$null_model
  dispersion <- numeric(length(x))
  dispersion[which(groups == 0)] <- theta[1]
  dispersion[which(groups == 1)] <- theta[2]
  if(null_model) {
    mu <- rep(exp(theta[3]), length(groups))
  } else {
    mu <- exp(theta[3] + groups*theta[4])
  }
  res = numeric(length(x))
  for(i in 1:length(x)) {
    log.gamma <- 0
    if(x[i] > 0) {
      log.gamma <- sum(sapply(1:x[i], function(x) log(x)))
    }
    res[i] = res[i] - log.gamma - (dispersion[i]^(-1)+x[i])*log(1+dispersion[i]*size_factor[i]*mu[i]) + x[i]*log(size_factor[i]*mu[i]) + x[i]*log(dispersion[i])
    if(x[i] > 0) {
      for(j in 0:(x[i]-1)) {
        res[i] = res[i] + log(j + dispersion[i]^(-1))
      }
    }
  }
  -sum(res)
}

#' Attempt to recreate the observed counts for a given gene from the fitted model (giving a means to eyeball
#' model goodness-of-fit)
#'
#' @param theta vector of parameters (see details)
#' @param other_params named list of counts, group assignments, size factor (per-sample total counts), and H0 flag
#' @details Parameter vector contains (in this order) control group dispersion, treatment group dispersion,
#' focal gene intercept, and (if not evaluating null model) focal gene treatment group effect.
#' @return a set of predicted counts from the fitted model
#' @export
predict.NBID <- function(theta, other_params) {
  x <- other_params$x
  groups <- other_params$group
  size_factor <- other_params$size_factor # this means mu is learned as a proportion
  #size_factor <- rep(mean(x), length(x)) # this means mu is learned as a scale
  null_model <- other_params$null_model
  dispersion <- groups
  dispersion[which(dispersion == 0)] <- theta[1]
  dispersion[which(dispersion == 1)] <- theta[2]
  if(null_model) {
    mu <- rep(exp(theta[3]), length(groups))
  } else {
    mu <- exp(theta[3] + groups*theta[4])
  }
  predicted <- sapply(1:length(x), function(i) {
    rnbinom(1, mu = mu[i]*size_factor[i], size = 1/dispersion[i])
  })
  return(predicted)
}

#' Calculate gradients for NBID model
#'
#' @param theta vector of parameters (see details)
#' @param other_params named list of counts, group assignments, size factor (per-sample total counts), and H0 flag
#' @details Parameter vector contains (in this order) control group dispersion, treatment group dispersion,
#' focal gene intercept, and (if not evaluating null model) focal gene treatment group effect.
#' @return gradient vector
#' @export
#' @examples
#' theta <- c(runif(2, min = 0.001, max = 1), runif(2, min = -1, max = 1))
#' other_params <- list(x = observed_counts, groups = sample_group_assignments, null_model = FALSE)
#' grad_vector <- gradient.NBID(theta, other_params)
gradient.NBID <- function(theta, other_params) {
  x <- other_params$x
  groups <- other_params$group
  size_factor <- other_params$size_factor # this means mu is learned as a proportion
  #size_factor <- rep(mean(x), length(x)) # this means mu is learned as a scale
  null_model <- other_params$null_model
  dispersion <- groups
  dispersion[which(dispersion == 0)] <- theta[1]
  dispersion[which(dispersion == 1)] <- theta[2]
  values <- numeric(length(theta))
  if(null_model) {
    values <- numeric(3)
    mu <- rep(exp(theta[3]), length(groups))
  } else {
    values <- numeric(4)
    mu <- exp(theta[3] + groups*theta[4])
  }
  # dispersion group 0, dispersion group 1, intercept, group 1 effect
  for(i in 1:length(x)) {
    a.i <- 1 + dispersion[i]*size_factor[i]*mu[i]
    a.inv <- 1/dispersion[i]
    a.inv.sq <- a.inv^2
    ln.a <- log(a.i)
    mu.diff <- x[i] - size_factor[i]*mu[i]
    b <- 0
    if(x[i] > 0) {
      for(j in 0:(x[i]-1)) {
        b <- b + 1/(j + a.inv)
      }
    }
    if(groups[i] == 0) {
      dispersion_update_idx <- 1
    } else {
      dispersion_update_idx <- 2
    }
    values[dispersion_update_idx] <- values[dispersion_update_idx] + a.inv.sq*(ln.a - b) + mu.diff/(dispersion[i]*a.i)
    values[3] <- values[3] + mu.diff/a.i
    if(!null_model) {
      values[4] <- values[4] + (groups[i]*mu.diff)/a.i
    }
  }
  -values
}

#' Optimize parameters for NBID model
#'
#' @param theta vector of parameters (see details)
#' @param other_params named list of counts, group assignments, and H0 flag
#' @details Parameter vector contains (in this order) control group dispersion, treatment group dispersion,
#' focal gene intercept, and (if not evaluating null model) focal gene treatment group effect.
#' @return vector of estimated parameter values
#' @import nloptr
#' @export
#' @examples
#' theta <- c(runif(2, min = 0.001, max = 1), runif(2, min = -1, max = 1))
#' other_params <- list(x = observed_counts, groups = sample_group_assignments, null_model = FALSE)
#' theta_hat <- optimize.NBID(theta, other_params)
optimize.NBID <- function(theta, other_params) {
  opts <- list(
    "algorithm" = "NLOPT_LD_LBFGS",
    "xtol_rel" = 1.0e-8,
    "maxeval" = 1e4)
  lb <- c(rep(1e-05, 2), rep(-1e05, length(theta)-2))
  ub <- c(rep(10, 2), rep(1e05, length(theta)-2))
  resultNLOPT <- nloptr::nloptr(theta,
                                eval_f = negativeLL.NBID,
                                eval_grad_f = gradient.NBID,
                                opts = opts,
                                lb = lb,
                                ub = ub,
                                other_params = other_params)
  resultNLOPT$solution
}

#' Calculate p-value from likelihood ratio test over NBID full and null models
#'
#' @param theta_H0 estimated null parameter set
#' @param theta_H1 estimated full parameter set
#' @param other_params named list of counts and group assignments
#' @return p-value from two-tailed chi-squared test
#' @export
#' @examples
#' other_params <- list(x = observed_counts, groups = sample_group_assignments)
#' theta_H0 <- c(runif(2, min = 0.001, max = 1), runif(1, min = -1, max = 1))
#' other_params_H0 <- other_params
#' other_params_H0$null_model <- TRUE
#' theta_H0_fitted <- optimize.NBID(theta_H0, other_params_H0)
#' theta_H1 <- c(theta_H0, runif(1, min = -1, max = 1))
#' other_params_H1 <- other_params
#' other_params_H1$null_model <- TRUE
#' theta_H1_fitted <- optimize.NBID(theta_H1, other_params_H1)
#' pval <- calc_LRT(theta_H0_fitted, theta_H1_fitted, other_params)
calc_LRT.NBID <- function(theta_H0, theta_H1, other_params) {
  other_params$null_model <- TRUE
  ll_H0 <- -negativeLL.NBID(theta_H0, other_params)
  other_params$null_model <- FALSE
  ll_H1 <- -negativeLL.NBID(theta_H1, other_params)
  lrt_val <- -2*(ll_H0 - ll_H1)
  pchisq(lrt_val, df = abs(length(theta_H0) - length(theta_H1)), lower.tail=FALSE)
}
