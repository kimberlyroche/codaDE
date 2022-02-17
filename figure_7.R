source("path_fix.R")

library(ggplot2)
library(MASS)
library(cowplot)

# Read in NB GLM calls on absolute abundances and DESeq2 calls on relative 
# abundances for Klein et al., which shows huge false positive rates for all
# methods.

dataset_name <- "Owens"
method <- "scran"
call <- "FP"

data <- readRDS(paste0("output/filtered_data_",dataset_name,"_threshold1.rds"))

calls_obj <- readRDS(paste0("output/real_data_calls/no_norm/calls_oracle_",method,"_", dataset_name, "_threshold1.rds"))

# Parse data using code in `validate.R` before proceeding...

fgroups <- factor(data$groups)
palette <- c("#f25c4b", "#3199e8")

wrong_feature <- sample(which(calls_obj$rates[[paste0(call,"_calls")]] == TRUE), size = 1)
A_idx <- sample(which(fgroups == levels(fgroups)[1]))
B_idx <- sample(which(fgroups == levels(fgroups)[2]))
plot_df <- data.frame(x = 1:length(A_idx),
                      y = data$absolute[A_idx,wrong_feature],
                      condition = levels(fgroups)[1])
plot_df <- rbind(plot_df,
                 data.frame(x = (length(A_idx)+1):(length(A_idx)+length(B_idx)),
                            y = data$absolute[B_idx,wrong_feature],
                            condition = levels(fgroups)[2]))
p1 <- ggplot(plot_df, aes(x = x, y = y, fill = condition)) +
  geom_point(size = 2, shape = 21) +
  scale_fill_manual(values = palette) +
  theme_bw() +
  labs(x = "sample index",
       y = "nominal abundance")

plot_df <- data.frame(x = 1:length(A_idx),
                      y = data$relative[A_idx,wrong_feature],
                      condition = levels(fgroups)[1])
plot_df <- rbind(plot_df,
                 data.frame(x = (length(A_idx)+1):(length(A_idx)+length(B_idx)),
                            y = data$relative[B_idx,wrong_feature],
                            condition = levels(fgroups)[2]))
p2 <- ggplot(plot_df, aes(x = x, y = y, fill = condition)) +
  geom_point(size = 2, shape = 21) +
  scale_fill_manual(values = palette) +
  theme_bw() +
  labs(x = "sample index",
       y = "relative abundance")

if(method == "DESeq2") {
  scaled_counts <- scaled_counts_DESeq2(data$relative, data$groups, 0.5)
} else if(method == "scran") {
  scaled_counts <- scaled_counts_scran(data$relative, data$groups, 0.5)
} else {
  scaled_counts <- scaled_counts_ALDEx2(data$relative, 0.5)
}

plot_df <- data.frame(x = 1:length(A_idx),
                      y = scaled_counts[A_idx,wrong_feature],
                      condition = levels(fgroups)[1])
plot_df <- rbind(plot_df,
                 data.frame(x = (length(A_idx)+1):(length(A_idx)+length(B_idx)),
                            y = scaled_counts[B_idx,wrong_feature],
                            condition = levels(fgroups)[2]))
p3 <- ggplot(plot_df, aes(x = x, y = y, fill = condition)) +
  geom_point(size = 2, shape = 21) +
  scale_fill_manual(values = palette) +
  theme_bw() +
  labs(x = "sample index",
       y = paste0(method, "-scaled abundance"))

p <- plot_grid(p1 + theme(legend.position = "none"),
               p2 + theme(legend.position = "none"),
               p3 + theme(legend.position = "none"), ncol = 3,
               labels = c("a", "b", "c"),
               label_size = 18,
               scale = 0.92)
p

ggsave(paste0("output/images/", dataset_name, "_",call,"_example_",method,".png"),
       p,
       dpi = 100,
       units = "in",
       height = 3,
       width = 9)

# Re-run the test to make sure this is "not differential" by the NB GLM
# call_DA_NB(data = data$absolute[,wrong_feature,drop=FALSE], fgroups)
