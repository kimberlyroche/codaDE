library(ggplot2)
library(gridExtra)
library(codaDE)
# library(wesanderson)

args = commandArgs(trailingOnly=TRUE)
label <- args[1]

data_dir <- file.path("simulated_analyses", label)
results_file <- "results.tsv"

data <- read.table(file.path(data_dir, results_file), header = TRUE)

# calculate error rates of interest
data$tpr <- data$tp / (data$tp + data$fn)
data$fpr <- data$fp / (data$fp + data$tn)
data$proportion_da <- round(data$prop_da, 2)
data <- data[data$prop_da %in% c(0.2, 0.4, 0.8),]
data <- data[data$sf_corr %in% c(0.1, 0.5, 0.9),]

data$prop_da <- as.factor(data$prop_da)
data$sf_corr <- as.factor(data$sf_corr)

# filter out 0/0 that occasionally crops up in filtered ALR w/ low number of features
# these are tp = 0, fn = 0
data <- data[!(is.nan(data$tpr)),]

for(ng in unique(data$p)) {
  subdata <- data[data$p == ng,]
  p1 <- ggplot(subdata, aes(x = fpr, y = tpr, color = prop_da)) +
    geom_point(size = 1) +
    scale_color_manual(values=c("#E3AF14", "#E36A14", "#E31B14")) +
    xlim(c(0, 1)) +
    ylim(c(0, 1)) +
    xlab("false positive rate") +
    ylab("true postitive rate") +
    ggtitle(paste0("p = ", ng)) +
    labs(color = "proportion genes\ndifferentially\nexpressed")
  ggsave(paste0("ROC_",ng,"_01.png"),
        p1, dpi = 150, units = "in", height = 4, width = 6)

  p2 <- ggplot(subdata, aes(x = fpr, y = tpr, color = sf_corr)) +
    geom_point(size = 1) +
    scale_color_manual(values=c("#11BF43", "#22a4d4", "#2225D4")) +
    xlim(c(0, 1)) +
    ylim(c(0, 1)) +
    xlab("false positive rate") +
    ylab("true postitive rate") +
    ggtitle(paste0("p = ", ng)) +
    labs(color = "size factor\ncorrelation")

  ggsave(paste0("ROC_",ng,"_02.png"),
          p2, dpi = 150, units = "in", height = 4, width = 6)
  #p <- grid.arrange(p1, p2, nrow = 1)
}
