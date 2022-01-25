source("path_fix.R")

library(ggplot2)
library(MASS)
library(cowplot)

# Read in NB GLM calls on absolute abundances and DESeq2 calls on relative 
# abundances for Klein et al., which shows huge false positive rates for all
# methods.

dataset_name <- "Hashimshony"

data <- readRDS(paste0("output/filtered_data_",dataset_name,"_threshold1.rds"))

calls_obj <- readRDS(paste0("output/real_data_calls/no_norm/calls_oracle_DESeq2_", dataset_name, "_threshold1.rds"))

# Parse data using code in `validate.R` before proceeding...

fgroups <- factor(data$groups)

wrong_feature <- sample(which(calls_obj$rates$FP_calls == TRUE), size = 1)
A_idx <- which(fgroups == levels(fgroups)[1])
B_idx <- which(fgroups == levels(fgroups)[2])
plot_df <- data.frame(x = 1:length(A_idx),
                      y = data$absolute[A_idx,wrong_feature],
                      condition = levels(fgroups)[1])
plot_df <- rbind(plot_df,
                 data.frame(x = (length(A_idx)+1):(length(A_idx)+length(B_idx)),
                            y = data$absolute[B_idx,wrong_feature],
                            condition = levels(fgroups)[2]))
p1 <- ggplot(plot_df, aes(x = x, y = y, fill = condition)) +
  geom_point(size = 2, shape = 21) +
  theme_bw() +
  labs(x = "sample index",
       y = "absolute abundance")

plot_df <- data.frame(x = 1:length(A_idx),
                      y = data$relative[A_idx,wrong_feature],
                      condition = levels(fgroups)[1])
plot_df <- rbind(plot_df,
                 data.frame(x = (length(A_idx)+1):(length(A_idx)+length(B_idx)),
                            y = data$relative[B_idx,wrong_feature],
                            condition = levels(fgroups)[2]))
p2 <- ggplot(plot_df, aes(x = x, y = y, fill = condition)) +
  geom_point(size = 2, shape = 21) +
  theme_bw() +
  labs(x = "sample index",
       y = "relative abundance")

scaled_counts <- scaled_counts_DESeq2(data$relative, data$groups, 0.5)

plot_df <- data.frame(x = 1:length(A_idx),
                      y = scaled_counts[A_idx,wrong_feature],
                      condition = levels(fgroups)[1])
plot_df <- rbind(plot_df,
                 data.frame(x = (length(A_idx)+1):(length(A_idx)+length(B_idx)),
                            y = scaled_counts[B_idx,wrong_feature],
                            condition = levels(fgroups)[2]))
p3 <- ggplot(plot_df, aes(x = x, y = y, fill = condition)) +
  geom_point(size = 2, shape = 21) +
  theme_bw() +
  labs(x = "sample index",
       y = "scaled abundance")

p <- plot_grid(p1, p2, p3, ncol = 3)
show(p)

# ggsave(paste0("output/images/", dataset_name, "_FP_example.png"),
#        p,
#        dpi = 100,
#        units = "in",
#        height = 4,
#        width = 5)

# Re-run the test to make sure this is "not differential" by the NB GLM
call_DA_NB(data = data$absolute[,wrong_feature,drop=FALSE], fgroups)
