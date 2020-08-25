library(ggplot2)
library(gridExtra)
library(codaDE)
library(wesanderson)

# options:
# (1) counts
# (2) logratios
measure <- "counts"
# options:
# (1) 0   -- no filtering of lowly abundant features
# (1) 1   -- omission of calls on features with mean abundance < 1
min_abundance <- 0
# options:
# (1) TRUE  -- use NB GLM + LRT for DE calling
# (2) FALSE -- use log-LM + permutations for DE calling
NB_DE <- FALSE

# output_file <- "results_singlecellRNAseq.tsv"
output_file <- "results_singlecellRNAseq_smallbatch.tsv"

data <- read.table(file.path("output", output_file), header = F)
colnames(data) <- c("run_label", "p", "prop_de", "sfcorr", "measure", "NB_DE",
                    "min_abundance", "tp", "fp", "tn", "fn", "feature_evaluated")
# subset to a condition of interest
data <- data[data$measure == measure & data$min_abundance == min_abundance & data$NB_DE == NB_DE,]

# calculate error rates of interest
data$tpr <- data$tp / (data$tp + data$fn)
data$fpr <- data$fp / (data$fp + data$tn)
data$prop_de <- round(data$prop_de, 2)
data$prop_de <- as.factor(data$prop_de)
data$sfcorr <- as.factor(data$sfcorr)

# filter out 0/0 that occasionally crops up in filtered ALR w/ low number of features
# these are tp = 0, fn = 0
data <- data[!(is.nan(data$tpr)),]

for(ng in unique(data$p)) {
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
  ggsave(paste0("ROC_",ng,"_",measure,"_min",min_abundance,"_01.png"),
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

  ggsave(paste0("ROC_",ng,"_",measure,"_min",min_abundance,"_02.png"),
          p2, dpi = 150, units = "in", height = 4, width = 6)
  #p <- grid.arrange(p1, p2, nrow = 1)

