library(ggplot2)

rm(list = ls())

cvs <- function(counts) {
  # Coefficient of variation
  sigmas <- apply(counts, 1, sd)
  mus <- apply(counts, 1, mean)
  retain_idx <- which(mus > 0)
  cvs <- sigmas[retain_idx] / mus[retain_idx]
  return(cvs)
}

# Plot binned coefficient of variation distribution for lowly, moderately and highly expressed genes.

sim_cartoon_foldchange <- FALSE
# data_type <- list("yeast", "C:/Users/kim/Documents/codaDE/data/Athanasiadou_2021/S1CodeandData/yeast_parsed.rds")
# data_type <- list("ciona", "C:/Users/kim/Documents/codaDE/data/Athanasiadou_2021/S1CodeandData/ciona_parsed.rds")
# data_type <- list("absolute1", "C:/Users/kim/Documents/codaDE/data/Barlow_2020/absolute_parsed.rds")
data_type <- list("absolute2", "C:/Users/kim/Documents/codaDE/data/Morton_2019/Github/absolute_parsed.rds")
data_obj <- readRDS(data_type[[2]])
counts <- data_obj$counts

if(data_type[[1]] == "ciona") {
  # Ciona counts are super low; stretch these (proportionally) into a more useful range
  counts <- counts * 100
}

conditions <- data_obj$groups
condition_A <- as.character(levels(conditions)[1])
if(length(levels(conditions)) == 2) {
  condition_B <- as.character(levels(conditions)[2])
} else {
  condition_B <- as.character(levels(conditions)[length(levels(conditions))])
}

# Plot overall distributions (condition A)
plot_data <- data.frame(log_counts = log(c(counts[,conditions == condition_A] + 1)))

ggplot(plot_data, aes(x = log_counts)) +
  geom_histogram(bins = 30) +
  xlim(c(-0.5, 15))
ggsave(paste0("density_",data_type[[1]],".png"), units = "in", dpi = 100, height = 5, width = 6)

if(sim_cartoon_foldchange) {
  # Replace observed counts in condition B with some random fold change (0.5x - 2x) to
  # see how this changes distributions
  for(i in 1:nrow(counts)) {
    counts[i,7:9] <- counts[i,1:3] * runif(n = 1, min = 0.5, max = 2)
  }
}

# Plot totals
plot_data <- data.frame(total = colSums(counts[,conditions == condition_A]), type = condition_A)
plot_data <- rbind(plot_data, data.frame(total = colSums(counts[,conditions == condition_B]), type = condition_B))
plot_data$sample_index <- 1:nrow(plot_data)
plot_data$type <- factor(plot_data$type, levels = c(condition_A, condition_B))

ggplot(plot_data, aes(x = sample_index, y = total, total, fill = type)) +
  geom_bar(stat = "identity")
ggsave(paste0("barplot_",data_type[[1]],".png"), units = "in", dpi = 100, height = 5, width = 6)

mean_expression_A <- rowMeans(counts[,conditions == condition_A])
mean_expression_B <- rowMeans(counts[,conditions == condition_B])

# binned_expression <- cut(mean_expression_A, breaks = c(0, 1, 10, 20, 100, 1000, max(mean_expression_A)))
# table(binned_expression)

breaks <- quantile(mean_expression_A, probs = c(0.5, 0.95))
breaks

low_genes <- mean_expression_A < breaks[1]
med_genes <- mean_expression_A >= breaks[1] & mean_expression_A < breaks[2]
high_genes <- mean_expression_A >= breaks[2]

cv_low <- cvs(counts[low_genes, conditions %in% c(condition_A, condition_B)])
cv_med <- cvs(counts[med_genes, conditions %in% c(condition_A, condition_B)])
cv_high <- cvs(counts[high_genes, conditions %in% c(condition_A, condition_B)])

plot_data <- data.frame(cv = cv_low, type = "low")
plot_data <- rbind(plot_data, data.frame(cv = cv_med, type = "med"))
plot_data <- rbind(plot_data, data.frame(cv = cv_high, type = "high"))
plot_data$type <- factor(plot_data$type, levels = c("low", "med", "high"))

ggplot(plot_data, aes(x = cv, color = type)) +
  geom_density(size = 2) +
  xlab("coefficient of variation")
if(sim_cartoon_foldchange) {
  ggsave(paste0("cv_",data_type[[1]],"_sim.png"), units = "in", dpi = 100, height = 5, width = 6)
} else {
  ggsave(paste0("cv_",data_type[[1]],".png"), units = "in", dpi = 100, height = 5, width = 6)
}


