library(MASS)

# Build baseline expression profile

# Function to estimate mean and dispersion of a count vector assuming an NB model
fit_gene <- function(idx, expression_mat, condition_columns) {
  expression <- unname(expression_mat[idx,condition_columns])
  if(var(expression) == 0) {
    expression <- expression + sample(c(1, rep(0, length(expression) - 1)))
  }
  fit <- tryCatch({
    suppressWarnings(glm.nb(y ~ 1, data = data.frame(y = expression)))
  },
  error = function(cond) { list() }
  )
  if(length(fit) > 0) {
    est_mu <- unname(exp(fit$coefficients[1]))
    est_size <- fit$theta
  } else {
    est_mu <- mean(expression)
    est_size <- 0
  }
  if(est_size < 1) {
    est_size <- 1 # What's a good default?
  }
  return(list(idx = idx, x = expression, mu = est_mu, size = est_size))
}

# data_type <- list("yeast", "C:/Users/kim/Documents/codaDE/data/Athanasiadou_2021/S1CodeandData/yeast_parsed.rds")
# data_type <- list("ciona", "C:/Users/kim/Documents/codaDE/data/Athanasiadou_2021/S1CodeandData/ciona_parsed.rds")
# data_type <- list("absolute1", "C:/Users/kim/Documents/codaDE/data/Barlow_2020/absolute_parsed.rds")
data_type <- list("absolute2", "C:/Users/kim/Documents/codaDE/data/Morton_2019/Github/absolute_parsed.rds")

data_obj <- readRDS(data_type[[2]])
counts <- data_obj$counts
if(data_type[[1]] == "ciona") {
  counts <- counts * 100
}

conditions <- data_obj$groups
condition_A <- as.character(levels(conditions)[1])
if(length(levels(conditions)) == 2) {
  condition_B <- as.character(levels(conditions)[2])
} else {
  condition_B <- as.character(levels(conditions)[length(levels(conditions))])
}

n_genes <- nrow(counts)
baseline_distribution <- matrix(NA, n_genes, 2)
for(i in 1:n_genes) {
  if(i %% 1000 == 0) {
    cat("Iteration",i,"...\n")
  }
  fit_obj <- fit_gene(i, counts, conditions == condition_A)
  baseline_distribution[i,] <- c(fit_obj$mu, fit_obj$size)
}

# n_genes <- 100 # testing
differential_distribution <- NULL
for(i in 1:n_genes) {
  if(i %% 1000 == 0) {
    cat("Iteration",i,"...\n")
  }
  expr1 <- unname(counts[i, conditions == condition_A])
  expr2 <- unname(counts[i, conditions == condition_B])
  if(var(expr1) == 0 & var(expr2) == 0) {
    next;
  }
  fit <- suppressWarnings(glm.nb(y ~ condition, data = data.frame(y = c(expr1, expr2),
                                           condition = as.factor(c(rep(condition_A, length(expr1)), rep(condition_B, length(expr2)))))))
  pval <- coef(summary(fit))[2,4]
  if(pval < 0.05) {
    fit_obj <- fit_gene(i, counts, conditions == condition_B)
    if(is.null(differential_distribution)) {
      differential_distribution <- c(baseline_distribution[i,], fit_obj$mu, fit_obj$size)
    } else {
      differential_distribution <- rbind(differential_distribution, c(baseline_distribution[i,], fit_obj$mu, fit_obj$size))
    }
  }
}
rownames(differential_distribution) <- NULL

cat("Proportion apparent DE:", round(nrow(differential_distribution)/n_genes, 2), "\n")

saveRDS(list(baseline_distribution = baseline_distribution, differential_distribution = differential_distribution),
        file = paste0("empirical_distro_",data_type[[1]],".rds"))

# Simulate differential expression by drawing a baseline mean abundance and matching that to a (random but similar)
# mean abundance in the differential expression distribution.

rm(list = ls())

data_obj <- readRDS("empirical_distro_absolute1.rds")

differential_baselines <- data_obj$differential_distribution[,1]
breaks <- quantile(differential_baselines, probs = seq(from = 0, to = 1, length.out = 5))
# Fix for R's cut function being kind of shitty
breaks[[1]] <- -Inf
breaks[[length(breaks)]] <- Inf

baseline_idx <- sample(1:nrow(data_obj$baseline_distribution), size = 1)
baseline_mu <- data_obj$baseline_distribution[baseline_idx,1]
baseline_size <- data_obj$baseline_distribution[baseline_idx,2]

binned_expression <- cut(c(baseline_mu, differential_baselines), breaks = breaks)
baseline_mu_binned <- binned_expression[1]
binned_expression <- binned_expression[2:length(binned_expression)]
similar_mus <- which(baseline_mu_binned == binned_expression)

differential_idx <- sample(similar_mus, size = 1)
differential_mu <- data_obj$differential_distribution[differential_idx,3]
differential_size <- data_obj$differential_distribution[differential_idx,4]

expr1 <- rnbinom(100, mu = baseline_mu, size = baseline_size)
expr2 <- rnbinom(100, mu = differential_mu, size = differential_size)
plot(c(expr1, expr2))









