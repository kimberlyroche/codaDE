library(dplyr)
library(ggplot2)

sim_dir <- "simulated_analyses"
files <- list.files(sim_dir)

df <- data.frame(FPR = c(), method = c(), prop_da = c(), sf_corr = c())
for(file in files) {
	data <- readRDS(file.path(sim_dir, file))
	df <- rbind(df, data.frame(p = data$p,
                             FPR = data$FPR,
                             method = data$method,
                             prop_da = data$proportion_da,
                             sf_corr = data$size_factor_correlation))
}

# Conditions:
# (1) Sweep over number of features (1K, 5K, 10K, 20K), fixing DE at 50%, lib. sz. correlation at 0
# (2) Sweep over number of DE (20%, 40%, 60%, 80%), fixing features at 20K, lib. sz. correlation at 0
# (3) Sweep over number of DE (20%, 40%, 60%, 80%), fixing features at 20K, lib. sz. correlation at 0.5

for(condition in c(1,3)) {
  if(condition == 1) {
    # 1: Effect of feature (gene) number
    sub_df <- df[(df$sf_corr == 0 & df$method == "edgeR"),]
    p <- ggplot(sub_df, aes(x = as.factor(prop_da), y = FPR)) +
      geom_boxplot() +
      ylim(0, 0.75)
  } else if(condition == 2) {
    # 2: Library size normalization
    sub_df <- df[(df$sf_corr == 0 & df$method == "edgeR"),]
    p <- ggplot(sub_df, aes(x = as.factor(prop_da), y = FPR)) +
      geom_boxplot() +
      ylim(0, 0.25)
    ggsave(paste0("test_",condition,".png"), p, units = "in", dpi = 100, height = 5, width = 5)
  } else {
    # 3: Faithful library size preservation
    sub_df <- df[(df$prop_da == 0.5 & df$method == "NB"),]
    p <- ggplot(sub_df, aes(x = as.factor(sf_corr), y = FPR)) +
      geom_boxplot() +
      ylim(0, 0.25)
    ggsave(paste0("test_",condition,".png"), p, units = "in", dpi = 100, height = 5, width = 5)
  }
}



# Calculate fold change in total abundance
sub_df <- df[(df$sf_corr == 0 & df$method == "edgeR"),]
m1 <- mean(data$library_size.abundances[1:100])
m2 <- mean(data$library_size.abundances[101:200])
cat("Approx. fold change in total abundance:",round(m2/m1, 2),"\n")
