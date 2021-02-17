library(ggplot2)
library(gridExtra)
library(codaDE)
library(RColorBrewer)
library(tidyverse)

generate_highcontrast_palette <- function(S) {
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  sample(getPalette(S))
}

n <- 20
p <- 10000
sim_data <- simulate_sequence_counts(n = n, p = p, k = 1, ref_data = "athanasiadou_ciona",
                                     sequencing_depth = 1e5, proportion_da = 0.5, library_size_correlation = 0,
                                     spike_in = FALSE, possible_fold_changes = NULL)

# Takes ~2min. with sequencing depth = 1e5 and p = 20K
calls <- call_DA_edgeR(sim_data, call_abundances = FALSE, normalization_method = NULL)
str(calls)

de <- which(calls < 0.05)
not_de <- setdiff(1:length(calls), de)

FP <- setdiff(de, sim_data$da_assignment)
FN <- setdiff(not_de, sim_data$da_assignment)

# FPR
length(FP) / ncol(sim_data$abundances)

# FNR
length(FN) / ncol(sim_data$abundances)

# Delta total abundance
# plot(rowSums(sim_data$abundances))
# plot(rowSums(sim_data$observed_counts))
# 
# idx <- sample(FP, size = 1)
# p1 <- ggplot(data.frame(y = c(sim_data$abundances[,idx]),
#                         x = 1:(n*2),
#                         tissue = as.factor(c(rep("A", n), rep("B", n)))), aes(x = x, y = y, color = tissue)) +
#   geom_point(size = 2) +
#   xlab("sample index") +
#   ylab("abundance")
# p2 <- ggplot(data.frame(y = c(sim_data$observed_counts[,idx]),
#                         x = 1:(n*2),
#                         tissue = as.factor(c(rep("A", n), rep("B", n)))), aes(x = x, y = y, color = tissue)) +
#   geom_point(size = 2) +
#   xlab("sample index") +
#   ylab("abundance")
# grid.arrange(grobs = list(p1, p2), ncol = 2)

# Plot differences in total abundance between conditions
# data <- data.frame(total_counts = rowSums(sim_data$abundances), condition = as.factor(sim_data$groups))
# 
# ggplot(data, aes(x = total_counts, color = condition)) +
#   geom_density()
# 
# mean(data[data$condition == "1",]$total_counts) / mean(data[data$condition == "0",]$total_counts)

# For microbial plots, visualize stacked bar plots
# Is this amount of compositional variation realistic?

# Plot stacked bars -- relative abundance

# Eliminate very low proportion taxa
if(p <= 1000) {
  # Convert to proportions
  props <- t(apply(sim_data$abundances, 1, function(x) x/sum(x)))
  vlow_abundance <- (apply(props, 2, mean) < 0.01)
  props2 <- props[,!vlow_abundance]
  props2 <- cbind(props2, rowSums(props[,vlow_abundance]))
  props <- props2
} else {
  # Plot a subset of the composition
  props <- t(apply(sim_data$abundances[,sample(1:ncol(sim_data$abundances), size = 1000)], 1, function(x) x/sum(x)))
}

palette <- generate_highcontrast_palette(ncol(props))

data <- pivot_longer(cbind(sample = 1:nrow(props), as.data.frame(props)),
                     cols = !sample,
                     names_to = "OTU",
                     values_to = "relative_abundance")
ggplot(data, aes(fill = OTU, y = relative_abundance, x = sample)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = palette) +
  theme(legend.position = "none")

# Plot stacked bars -- absolute abundance

counts <- sim_data$abundances # samples x taxa
counts2 <- counts[,!vlow_abundance]
counts2 <- cbind(counts2, rowSums(counts[,vlow_abundance]))
counts <- counts2

data <- pivot_longer(cbind(sample = 1:nrow(counts), as.data.frame(counts)),
                     cols = !sample,
                     names_to = "OTU",
                     values_to = "abundance")

palette <- generate_highcontrast_palette(length(unique(data$OTU)))

ggplot(data, aes(fill = OTU, y = abundance, x = sample)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = palette) +
  theme(legend.position = "none")


