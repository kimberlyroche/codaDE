library(ggplot2)
library(gridExtra)
library(codaDE)
library(wesanderson)

measure <- "counts"
# options:
# (1) 0   -- no filtering of lowly abundant features
# (1) 1   -- omission of calls on features with mean abundance < 1
min_abundance <- 0

output_file <- "results_filter0.tsv"

data <- read.table(file.path("simulated_data", output_file), header = T)
# subset to a condition of interest
data <- data[data$filter_threshold == min_abundance,]

# length(data$filename %in% list.files(path = "simulated_data", pattern = "*.rds"))
# length(list.files(path = "simulated_data", pattern = "*.rds") %in% data$filename)

# calculate error rates of interest
data$tpr <- data$tp / (data$tp + data$fn)
data$fpr <- data$fp / (data$fp + data$tn)
data$proportion_da <- round(data$proportion_da, 2)
data <- data[data$proportion_da %in% c(0.2, 0.4, 0.8),]
data <- data[data$size_factor_correlation %in% c(0.1, 0.5, 0.9),]

data$proportion_da <- as.factor(data$proportion_da)
data$size_factor_correlation <- as.factor(data$size_factor_correlation)

# filter out 0/0 that occasionally crops up in filtered ALR w/ low number of features
# these are tp = 0, fn = 0
data <- data[!(is.nan(data$tpr)),]

for(ng in unique(data$p)) {
  subdata <- data[data$p == ng,]
  p1 <- ggplot(subdata, aes(x = fpr, y = tpr, color = proportion_da)) +
    geom_point(size = 1) +
    scale_color_manual(values=c("#E3AF14", "#E36A14", "#E31B14")) +
    xlim(c(0, 1)) +
    ylim(c(0, 1)) +
    xlab("false positive rate") +
    ylab("true postitive rate") +
    ggtitle(paste0("p = ", ng)) +
    labs(color = "proportion genes\ndifferentially\nexpressed")
  ggsave(paste0("ROC_",ng,"_min",min_abundance,"_01.png"),
        p1, dpi = 150, units = "in", height = 4, width = 6)

  p2 <- ggplot(subdata, aes(x = fpr, y = tpr, color = size_factor_correlation)) +
    geom_point(size = 1) +
    scale_color_manual(values=c("#11BF43", "#22a4d4", "#2225D4")) +
    xlim(c(0, 1)) +
    ylim(c(0, 1)) +
    xlab("false positive rate") +
    ylab("true postitive rate") +
    ggtitle(paste0("p = ", ng)) +
    labs(color = "size factor\ncorrelation")

  ggsave(paste0("ROC_",ng,"_min",min_abundance,"_02.png"),
          p2, dpi = 150, units = "in", height = 4, width = 6)
  #p <- grid.arrange(p1, p2, nrow = 1)
}

