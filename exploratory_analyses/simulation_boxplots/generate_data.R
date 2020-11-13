library(codaDE)
library(optparse)

option_list = list(
  make_option(c("--cond"), type = "numeric", default = NULL,
              help = "simulated condition ID", metavar = "numeric")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt$cond)) {
  stop("Missing simulated condition ID!")
}

if(opt$cond < 1 | opt$cond > 2) {
  stop("Invalid simulation condition!")
}

# Generate data for boxplots of simulated TPR and FPR.

# We want to simulate 2 conditions:
#  (1) Fixed library size correlation, increasing DE
#  (2) Fixed DE, decreasing library size correlation

if(opt$cond == 1) {
  prop_de_vec <- c(0.2, 0.4, 0.6)
  sf_corr_vec <- c(0.5)
} else {
  prop_de_vec <- c(0.4)
  sf_corr_vec <- c(0, 0.5)
}

# First, estimate an empirical fold change profile between tissues.
# Parse the GTEx data set.
GTEx <- readRDS(file.path("data", "GTEx_data", "parsed_GTEx.rds"))
GTEx_annot <- read.table(file.path("data", "GTEx_data", "GTEx_annotations.txt"), header = TRUE, sep = "\t")
samples_by_tissue <- list()
for(tissue in unique(GTEx_annot$SMTSD)) {
  samples_by_tissue[[tissue]] <- GTEx_annot$SAMPID[GTEx_annot$SMTSD == tissue]
}

# Pull a random pair of tissues.
tissue_pair <- sample(names(samples_by_tissue), size = 2, replace = FALSE)
cat("Using tissues:",tissue_pair[1],"x",tissue_pair[2],"\n")
# Calculate the mean expression in each tissue across all genes.
condition_means <- cbind(rowMeans(GTEx[,colnames(GTEx) %in% samples_by_tissue[[tissue_pair[1]]]]),
                         rowMeans(GTEx[,colnames(GTEx) %in% samples_by_tissue[[tissue_pair[2]]]]))
# Consider only genes that have substantially non-zero expression in at least one condition.
expr_lower_cutoff <- 1
keep <- which(condition_means[,1] >= expr_lower_cutoff | condition_means[,2] >= expr_lower_cutoff)
condition_means <- condition_means[keep,]

# Calculate and visualize fold change between non-zero genes. Log fold change closely resembles a Laplace distribution.
empirical_fold_change <- (condition_means[,1] + 0.001) / (condition_means[,2] + 0.001) # tiny addend prevents infinities

# Require a minimum fold change of 2 (or 1/2) for differentially expressed genes.
emp_fc <- empirical_fold_change[empirical_fold_change >= 2 | empirical_fold_change <= 0.5]

# This takes ~ 7 seconds for p = 1000 x n = 100.
# The time increase seems to be linear such that I'd expect about 5 min. for p = 20000 x n = 200.
# Note: We only achieve about 1/2 to 2/3 the differential expression we request. This is because even a 2-fold change
#       can be imperceptible from noise.
results <- data.frame(prop_de = c(), sf_corr = c(), TPR = c(), FPR = c(), prop_de_detectable = c(), runtime = c())
p <- 20000
n <- 200
for(prop_de in prop_de_vec) {
  for(sf_corr in sf_corr_vec) {
    cat("Evaluating % DE",round(prop_de*100),"x library size correlation",round(sf_corr,2),"\n")
    # fc <- c(50)
    fc <- emp_fc
    data <- simulate_singlecell_RNAseq(p = p, n = n, k = 1, proportion_da = prop_de, size_factor_correlation = sf_corr,
                                       spike_in = FALSE, possible_fold_changes = fc)
    ptm <- proc.time()
    res <- evaluate_DA(data, alpha = 0.05, use_ALR = FALSE, filter_abundance = 0, method = "edgeR", save_output = FALSE)
    runtime <- proc.time() - ptm
    runtime <- unname(runtime[3])
    results <- rbind(results,
                     data.frame(prop_de = prop_de, sf_corr = sf_corr, TPR = res$TPR, FPR = res$FPR,
                                prop_de_detectable = res$no_features_detectable/p, runtime = runtime))
  }
}
saveRDS(results, file.path("output", paste0("results_",opt$cond,"_",as.numeric(Sys.time()),".rds")))
