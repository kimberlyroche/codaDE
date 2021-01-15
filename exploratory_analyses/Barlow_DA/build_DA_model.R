# Generate empirical sample from Barlow et al. (2020) data -- as GTEx

library(MASS)
library(ggplot2)

# Copied from build_DE_model_1.R

fit_nb_model <- function(x) {
  fit <- tryCatch({
    glm.nb(y ~ 1, data.frame(y = x))
  },
  error = function(cond) { list() },
  warning = function(cond) { list() },
  finally = {}
  )
  if(length(fit) > 0) {
    est_mu <- unname(exp(fit$coefficients[1]))
    est_size <- fit$theta
    return(c(est_mu, est_size))
  } else {
    est_mu <- mean(x)
    # Estimate the variance this way
    var_x <- var(x)
    if(var_x == 0) {
      var_x <- 0.1
    }
    est_size <- (est_mu^2) / (var_x - est_mu)
    return(c(est_mu, est_size))
  }
}

est_mu_size <- function(x) {
  if(sum(x) == 0) {
    return(c(0.5, 1))
  } else {
    params <- fit_nb_model(x)
    if(params[1] < 1) {
      params[1] <- 0.5
    }
    # if(is.na(params[2])) {
    #   params[2] <- 1
    # }
    return(params)
  }
}

base.path <- "C:/Users/kim/Desktop/Barlow_2020"

data <- read.table(file.path(base.path, "Absolute_Abundance_Table.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Clean up data
md_labels <- c("Diet", "Site", "Day", "mouse", "Cage")
metadata <- data[,md_labels]
counts <- data[,!(colnames(data) %in% md_labels)]
counts <- counts[,2:ncol(counts)]
tax <- colnames(counts)
counts <- as.matrix(counts)
colnames(counts) <- NULL
rownames(counts) <- NULL

# Counts is initially 103 samples x 142 taxa

# Pull stool samples from days 4, 7, 10 into
keto_idx <- which(metadata$Diet == "Keto" & metadata$Day %in% c(4,7,10) & metadata$Site == "Stool")
ctrl_idx <- which(metadata$Diet == "Control" & metadata$Day %in% c(4,7,10) & metadata$Site == "Stool")

counts <- rbind(counts[ctrl_idx,], counts[keto_idx,])
labels <- c(rep("control", length(ctrl_idx)), rep("keto", length(keto_idx)))

# Eliminate all-zero taxa
absent_tax_idx <- which(colSums(counts) == 0)
counts <- counts[,-absent_tax_idx]
tax <- tax[-absent_tax_idx]

# Scale down by minimum observed abundance (?)
min_observed <- min(counts[which(counts > 0, arr.ind = TRUE)])
counts <- counts / min_observed

empirical_mu_size <- matrix(NA, ncol(counts), 4)
for(i in 1:ncol(counts)) {
  # i <- sample(1:ncol(counts), size = 1)
  p1 <- est_mu_size(round(counts[labels == "control",i]))
  p2 <- est_mu_size(round(counts[labels == "keto",i]))
  # Visualize
  # par(mfrow = c(1, 2))
  # ref_data <- c(counts[labels == "control",idx], counts[labels == "keto",idx])
  # sim_data <- c(rnbinom(length(ctrl_idx), mu = p1[1], size = p1[2]), rnbinom(length(keto_idx), mu = p2[1], size = p2[2]))
  # plot(ref_data, ylim = c(0, max(c(ref_data, sim_data))))
  # plot(sim_data, ylim = c(0, max(c(ref_data, sim_data))))
  empirical_mu_size[i,1:2] <- p1
  empirical_mu_size[i,3:4] <- p2
}

empirical_mu_size

# Estimate the delta total counts between conditions

totals_control <- rowSums(counts[labels == "control",])
totals_keto <- rowSums(counts[labels == "keto",])

data <- data.frame(x = c(totals_control, totals_keto), label = labels)
ggplot(data, aes(x = x, color = label)) +
  geom_density()

