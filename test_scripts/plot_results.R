library(ggplot2)

data <- read.table("output/results_RNAseq.tsv", header = F)
colnames(data) <- c("run_label", "p", "prop_de", "sfcorr", "tp", "fp", "tn", "fn")
data$tpr <- data$tp / (data$tp + data$fn)
data$fpr <- data$fp / (data$fp + data$tn)
data$prop_de <- as.factor(data$prop_de)
data$sfcorr <- as.factor(data$sfcorr)

p <- ggplot(data, aes(x = fpr, y = tpr, color = prop_de)) +
    geom_point(size = 1) +
    xlim(c(0, 1)) +
    ylim(c(0, 1)) +
    xlab("false positive rate") +
    ylab("true postitive rate") +
    labs(color = "proportion genes\ndifferentially\nexpressed") 
ggsave("RNAseq-like_color-DE.png", p, dpi = 150, units = "in", height = 8, width = 9)

p <- ggplot(data, aes(x = fpr, y = tpr, color = sfcorr)) +
    geom_point(size = 1) +
    xlim(c(0, 1)) +
    ylim(c(0, 1)) +
    xlab("false positive rate") +
    ylab("true postitive rate") +
    labs(color = "size factor\ncorrelation") 
ggsave("RNAseq-like_color-sizefactor.png", p, dpi = 150, units = "in", height = 8, width = 9)

