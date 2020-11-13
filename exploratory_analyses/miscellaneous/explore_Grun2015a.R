library(ggplot2)
library(gridExtra)

umi_data <- readRDS("../data/counts_conquerDB_Grun2015a.rds") # CEL-seq2 (?)

# percent zeros in this dataset
sum(umi_data$counts == 0)/(nrow(umi_data$counts)*ncol(umi_data$counts))

# the idea below is to recreate Fig. 1b and 1d from Townes et al.
# (1) a version of Fig. 1d: percent zeros x log mean abundance
# the authors downsampled to remove some variation due to total reads and (more importantly)
# were looking at *replicates*; I'm just interested in a high-level relationship (or lack)
# between mean log abundance and sparsity
mean_abundances <- rowMeans(umi_data$counts)
zero_percent <- apply(umi_data$counts, 1, function(x) {
  sum(x == 0)/length(x)
})
plot_df <- data.frame(x = log(mean_abundances + 0.5), y = zero_percent)
ggplot(plot_df, aes(x = x, y = y)) +
  geom_point()

# (2) this effectively repeats the percent zeros at the *sample* or cell level
total_counts <- colSums(umi_data$counts)
zero_percent_cell <- apply(umi_data$counts, 2, function(x) {
  sum(x == 0)/length(x)
})
plot_df <- data.frame(x = total_counts, y = zero_percent_cell)
ggplot(plot_df, aes(x = x, y = y)) +
  geom_point() +
  scale_x_continuous(trans='log10')

# tl;dr - zeros-per-cell are pretty strongly predicted by read depth
#         zeros-per-gene are pretty well predicted by average (log) expression level

# is there evidence of multiple modes in genes -- interpretable (maybe) as active/inactive transcription?
plot_df <- data.frame(x = log(mean_abundances + 0.5))
ggplot(plot_df, aes(x = x)) +
  geom_density()

# candidates <- which(mean_abundances > 5)
# g_idx <- sample(candidates)[1]
# print(unlist(g_idx))

# these genes are interesting: they illustrate that this data is obviously a mixture of cell types
genes <- c(4850, 11246)
for(g_idx in genes) {
  plot_df <- data.frame(x = log(umi_data$counts[g_idx,] + 0.5))
  p1 <- ggplot(plot_df, aes(x = x)) +
    geom_histogram() +
    xlab("log abundance") +
    ylab("counts")
  plot_df <- data.frame(x = 1:ncol(umi_data$counts), y = umi_data$counts[g_idx,])
  p2 <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_point() +
    xlab("sample index") +
    ylab("abundance")
  grid.arrange(p1, p2, ncol = 2)
}

mean_abundances <- rowMeans(umi_data$counts)
plot(density(log(mean_abundances + 0.5)))
max(mean_abundances)
