library(codaDE)
library(ggplot2)

plot_abundances <- function(idx, data, type = "abundances") {
  count_data <- data$abundances
  if(type == "observations") {
    count_data <- data$observed_counts
  }
  df <- data.frame(sample_id = 1:nrow(count_data),
                   count = count_data[,idx],
                   label = as.factor(data$groups))
  ggplot(df, aes(x = sample_id, y = count, color = label)) +
    geom_point()
}

plot_totals <- function(data, type = "abundances") {
  count_data <- data$abundances
  if(type == "observations") {
    count_data <- data$observed_counts
  }
  df <- data.frame(total = rowSums(count_data[1:n,]), label = "baseline")
  df <- rbind(df, data.frame(total = rowSums(count_data[(n+1):(2*n),]), label = "differential"))
  ggplot(df, aes(x = total, color = label)) +
    geom_density(size = 1)
}

n <- 250
prop_da <- 0.6
k <- 1

# data <- simulate_singlecell_RNAseq(n = n, k = k, proportion_da = prop_da, library_size_correlation = 0, possible_fold_changes = NULL)
data <- simulate_singlecell_RNAseq(n = n, k = k, proportion_da = prop_da, library_size_correlation = 0, possible_fold_changes = c(2,4,6,8))

avg_expr <- colMeans(data$abundances[1:n,data$da_assignment])
large_mean_idx <- data$da_assignment[sample(which(avg_expr > 2000), size = 1)]
plot_abundances(large_mean_idx, data)
plot_abundances(large_mean_idx, data, type = "observations")

small_mean_idx <- data$da_assignment[sample(which(avg_expr < 5), size = 1)]
plot_abundances(small_mean_idx, data)
plot_abundances(small_mean_idx, data, type = "observations")

plot_totals(data, type = "abundances")

