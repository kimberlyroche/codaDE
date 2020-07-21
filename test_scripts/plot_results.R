library(ggplot2)
library(gridExtra)
library(codaDE)
library(wesanderson)

plot_ROC <- TRUE
plot_DE <- FALSE

# options:
# (1) counts
# (2) logratios
measure <- "counts"
# options:
# (1) 0   -- no filtering of lowly abundant features
# (1) 1   -- omission of calls on features with mean abundance < 1
min_abundance <- 1
# data_type:
# (1) bulkRNAseq
# (2) singlecell_Haber (etc.)
data_type <- "bulkRNAseq"
no_genes <- c(100, 1000, 20000, 50000)

data_type <- "singlecell_Haber"
no_genes <- c(50000)

output_file <- paste0("results_", data_type, ".tsv")

# if(measure == "counts") {
#   output_file <- paste0(output_file, "RNAseq_")
# } else {
#   output_file <- paste0(output_file, "ALR_")
# }
# if(min_abundance == 0) {
#   output_file <- paste0(output_file, "nofilter.tsv")
# } else {
#   output_file <- paste0(output_file, "filter.tsv")
# }
# output files are:
# (1) results_RNAseq_nofilter.tsv   -- counts
# (1) results_RNAseq_filter.tsv     -- counts omitting tiny counts
# (2) results_ALR_nofilter.tsv      -- logratios
# (3) results_ALR_filter.tsv        -- logratios omitting tiny counts

if(plot_ROC) {
  data <- read.table(file.path("output", output_file), header = F)
  colnames(data) <- c("run_label", "p", "prop_de", "sfcorr", "measure", "min_abundance",
                      "tp", "fp", "tn", "fn", "feature_evaluated")
  # subset to a condition of interest
  data <- data[data$measure == measure & data$min_abundance == min_abundance,]
  data$tpr <- data$tp / (data$tp + data$fn)
  data$fpr <- data$fp / (data$fp + data$tn)
  # data$tnr <- data$tn / (data$tn + data$fp)
  data$prop_de <- round(data$prop_de, 2)
  data$prop_de <- as.factor(data$prop_de)
  data$sfcorr <- as.factor(data$sfcorr)

  # filter out 0/0 that occasionally crops up in filtered ALR w/ low number of features
  # these are tp = 0, fn = 0
  data <- data[!(is.nan(data$tpr)),]

  for(ng in no_genes) {
    subdata <- data[data$p == ng,]
    p1 <- ggplot(subdata, aes(x = fpr, y = tpr, color = prop_de)) +
      geom_point(size = 1) +
      scale_color_manual(values=c("#E3AF14", "#E36A14", "#E31B14")) +
      xlim(c(0, 1)) +
      ylim(c(0, 1)) +
      xlab("false positive rate") +
      ylab("true postitive rate") +
      ggtitle(paste0("p = ", ng)) +
      labs(color = "proportion genes\ndifferentially\nexpressed")
    ggsave(paste0(data_type,"_",ng,"_",measure,"_min",min_abundance,"_01.png"),
          p1, dpi = 150, units = "in", height = 4, width = 6)

    p2 <- ggplot(subdata, aes(x = fpr, y = tpr, color = sfcorr)) +
      geom_point(size = 1) +
      scale_color_manual(values=c("#11BF43", "#22a4d4", "#2225D4")) +
      xlim(c(0, 1)) +
      ylim(c(0, 1)) +
      xlab("false positive rate") +
      ylab("true postitive rate") +
      ggtitle(paste0("p = ", ng)) +
      labs(color = "size factor\ncorrelation")

    ggsave(paste0(data_type,"_",ng,"_",measure,"_min",min_abundance,"_02.png"),
           p2, dpi = 150, units = "in", height = 4, width = 6)
    #p <- grid.arrange(p1, p2, nrow = 1)
  }
}

# sanity check calls about differences
if(plot_DE) {
  p <- 20000
  data <- simulate_RNAseq(p = p, n = 250, proportion_de = (2/3), size_factor_correlation = 0)
  pvals.original <- c()
  pvals.resampled <- c()
  # the below takes about 1 min to run
  for(i in 1:(p/10)) {
    pvals.original <- c(pvals.original, call_DA(data, i, TRUE))
    pvals.resampled <- c(pvals.resampled, call_DA(data, i, FALSE))
  }
  # let's look at a false positive
  fp <- which((pvals.original > 0.05/p) & (pvals.resampled <= 0.05/p))
  fp_single <- sample(fp)[1]
  png("nonDE_sample.png", height = 400, width = 800)
  par(mfrow = c(1,2))
  plot(data$abundances[,fp_single], ylab = "original abundances")
  plot(data$observed_counts[,fp_single], ylab = "observed abundances")
  dev.off()
  cat("P-value abundance:", pvals.original[fp_single], "\n")
  cat("P-value observed counts:", pvals.resampled[fp_single], "\n")
}

