library(ggplot2)
library(gridExtra)
library(codaDE)

plot_ROC <- FALSE
plot_DE <- TRUE

if(plot_ROCE) {
  data <- read.table("output/results_RNAseq.tsv", header = F)
  colnames(data) <- c("run_label", "p", "prop_de", "sfcorr", "tp", "fp", "tn", "fn")
  data$tpr <- data$tp / (data$tp + data$fn)
  data$fpr <- data$fp / (data$fp + data$tn)
  data$prop_de <- as.factor(data$prop_de)
  data$sfcorr <- as.factor(data$sfcorr)

  no_genes <- c(100, 500, 20000)

  for(ng in no_genes) {
    subdata <- data[data$p == ng,]
    p1 <- ggplot(subdata, aes(x = fpr, y = tpr, color = prop_de)) +
      geom_point(size = 0.5) +
      xlim(c(0, 1)) +
      ylim(c(0, 1)) +
      xlab("false positive rate") +
      ylab("true postitive rate") +
      labs(color = "proportion genes\ndifferentially\nexpressed")

    p2 <- ggplot(subdata, aes(x = fpr, y = tpr, color = sfcorr)) +
      geom_point(size = 0.5) +
      xlim(c(0, 1)) +
      ylim(c(0, 1)) +
      xlab("false positive rate") +
      ylab("true postitive rate") +
      labs(color = "size factor\ncorrelation")

    p <- grid.arrange(p1, p2, nrow = 1)
    ggsave(paste0("RNAseq-like_p",ng,".png"), p, dpi = 150, units = "in", height = 4, width = 11)
  }
}

if(plot_DE) {
  p <- 500
  data <- simulate_RNAseq(p = p, n = 200, proportion_de = 0.2, size_factor_correlation = 0)
  for(i in sample(data$de_genes)[1]) {
    png("DE_sample.png", height = 400, width = 800)
    par(mfrow = c(1,2))
    plot(data$abundances[,i], ylab = "original abundances")
    plot(data$observed_counts[,i], ylab = "observed abundances")
    dev.off()
    cat("P-value abundance:", round(call_DE(data, i, TRUE), 3), "\n")
    cat("P-value observed counts:", round(call_DE(data, i, FALSE), 3), "\n")
  }
  for(i in sample(setdiff(1:p, data$de_genes))[1]) {
    png("nonDE_sample.png", height = 400, width = 800)
    par(mfrow = c(1,2))
    plot(data$abundances[,i], ylab = "original abundances")
    plot(data$observed_counts[,i], ylab = "observed abundances")
    dev.off()
    cat("P-value abundance:", round(call_DE(data, i, TRUE), 3), "\n")
    cat("P-value observed counts:", round(call_DE(data, i, FALSE), 3), "\n")
  }
}

