library(codaDE)
library(MASS)

# simulate some data
n <- 500
theta <- c(1, 1, 1, 0.5)
size_factors <- rep(100, n*2)
other_params <- list()
other_params$groups <- c(rep(0, n), rep(1, n))
other_params$x <- rnbinom(n*2, mu = size_factors*exp(theta[3] + other_params$groups*theta[4]), size = rep(1/theta[1], n*2)) # ignore separate group dispersion
other_params$size_factors <- rep(mean(other_params$x), length(other_params$x))
other_params$null_model <- FALSE

plot(other_params$x)

# compare negative log likelihood (R vs. C++)
negativeLL.NBID(theta, other_params)
negativeLL_NBID(theta, other_params$x, other_params$groups, other_params$size_factors, other_params$null_model)

# compare gradient given theta (R vs. C++)
gradient.NBID(theta, other_params)
gradient_NBID(theta, other_params$x, other_params$groups, other_params$size_factors, other_params$null_model)

theta_init <- runif(4, min = 0.1, max = 5)
optimize.NBID(theta_init, other_params)
optimize_NBID(theta_init, other_params$x, other_params$groups, other_params$size_factors, other_params$null_model)

data <- data.frame(y = other_params$x, group = as.factor(other_params$groups))
summary(m1 <- glm.nb(y ~ group, data = data))
