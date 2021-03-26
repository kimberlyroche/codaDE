library(codaDE)
library(tidyverse)

p <- 5000

# sim data
build_simulated_reference(p = p, log_mean = -1, log_var = 2, log_noise_var = 1)
ref_file <- file.path("data", paste0("DE_reference_simulated.rds"))
data_obj <- readRDS(ref_file)

# eyeball what observed log-counts would look like
ggplot(data.frame(x = log(round(data_obj$cond1) + 1)), aes(x = x)) +
  geom_histogram()

expr1 <- log(data_obj$cond1 + 1)
expr2 <- log(data_obj$cond2 + 1)
subset_idx <- sample(1:p, size = 100)
plot_bipartite_graph(expr1[subset_idx], expr2[subset_idx], alpha = 0.5)
