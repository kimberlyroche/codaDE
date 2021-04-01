# Run the internals of build_Morton_reference before running this

# Morton
labels <- subjects

# Athanasiadou
labels <- groups

pseudocount <- 1

# Bipartite graph, pairing samples for individual
expr1 <- c()
expr2 <- c()
for(i in 1:200) {
  # Pick a random subject
  subj <- sample(labels, size = 1)
  subj_idx <- which(labels == subj)
  # Pick a random feature
  feature <- sample(1:nrow(counts), size = 1)
  expr1 <- c(expr1, unname(unlist(counts[feature,subj_idx[1]])))
  expr2 <- c(expr2, unname(unlist(counts[feature,subj_idx[2]])))
}

plot_bipartite_graph(log(expr1 + pseudocount), log(expr2 + pseudocount))

# Bipartite graph, one individual; y-axis scaled to close the gap between observed and unobserved
subj <- sample(labels, size = 1)
subj_idx <- which(labels == subj)

subject_features <- counts[,c(subj_idx[1], subj_idx[2])]
subject_features <- subject_features[rowSums(subject_features)> 0,]

min_observed <- min(subject_features[which(subject_features > 0, arr.ind = TRUE)])
subject_features <- subject_features / min_observed

subset_n <- 100
subset_idx <- sample(1:nrow(subject_features), size = subset_n)
expr1 <- subject_features[subset_idx,1]
expr2 <- subject_features[subset_idx,2]

# plot_bipartite_graph(log(expr1 + pseudocount), log(expr2 + pseudocount))

# Stacked bar plots probably better exemplify the real story
data <- data.frame(abundance = c(expr1, expr2),
                   feature = rep(1:length(expr1), 2),
                   sample = c(rep("A", length(expr1)), rep("B", length(expr2))))
data$feature <- as.factor(data$feature)
data$sample <- as.factor(data$sample)
plot_stacked_bars(data)
