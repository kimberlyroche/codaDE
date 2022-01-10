source("path_fix.R")

library(tidyverse)
library(codaDE)
library(optparse)

option_list = list(
  make_option(c("--input"), type = "character", default = NULL,
              help = "input filename", metavar = "character"),
  make_option(c("--output"), type = "character", default = NULL,
              help = "output filename", metavar = "character"),
  make_option(c("--start"), type = "numeric", default = 1,
              help = "start row index", metavar = "numeric"),
  make_option(c("--end"), type = "numeric", default = Inf,
              help = "end row index", metavar = "numeric")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# ------------------------------------------------------------------------------
#  Validate input
# ------------------------------------------------------------------------------

input <- opt$input
output <- opt$output
start <- opt$start
end <- opt$end

if(!file.exists(input)) {
  stop("Input file not found!")
}
if(start < 1 | start > end) {
  stop("Invalid start row!")
}
if(end < 1) {
  stop("Invalid end row!")
}

output_fn <- file.path("temp", paste0(output, "_", start, "-", end, ".txt"))

# Initialize output
output_file <- file(output_fn)
writeLines(paste0(c("ID", "uuid", "partial_info", "type", "med_abs", "med_rel",
                    "TOTALS_C_FC", "TOTALS_C_D", "TOTALS_C_MAX_D",
                    "TOTALS_C_MED_D", "TOTALS_C_SD_D", "CORR_RA_MED",
                    "CORR_RA_SD", "CORR_RA_SKEW", "CORR_LOG_MED",
                    "CORR_LOG_SD", "CORR_LOG_SKEW", "CORR_CLR_MED",
                    "CORR_CLR_SD", "CORR_CLR_SKEW", "COMP_C_P0_A",
                    "COMP_C_P0_B", "COMP_C_P1_A", "COMP_C_P1_B",
                    "COMP_C_P5_A", "COMP_C_P5_B", "COMP_RA_P01_A",
                    "COMP_RA_P01_B", "COMP_RA_P1_A", "COMP_RA_P1_B",
                    "COMP_RA_P5_A", "COMP_RA_P5_B", "COMP_RA_MAX_A",
                    "COMP_RA_MED_A", "COMP_RA_SD_A", "COMP_RA_SKEW_A",
                    "COMP_RA_MAX_B", "COMP_RA_MED_B", "COMP_RA_SD_B",
                    "COMP_RA_SKEW_B", "COMP_C_ENT_A", "COMP_C_ENT_B",
                    "FW_RA_MAX_D", "FW_RA_MED_D", "FW_RA_SD_D",
                    "FW_RA_PPOS_D", "FW_RA_PNEG_D", "FW_RA_PFC05_D",
                    "FW_RA_PFC1_D", "FW_RA_PFC2_D", "FW_LOG_MAX_D",
                    "FW_LOG_MED_D", "FW_LOG_SD_D", "FW_LOG_PPOS_D",
                    "FW_LOG_PNEG_D", "FW_LOG_PFC05_D", "FW_LOG_PFC1_D",
                    "FW_LOG_PFC2_D", "FW_CLR_MAX_D", "FW_CLR_MED_D",
                    "FW_CLR_SD_D", "FW_CLR_PPOS_D", "FW_CLR_PNEG_D",
                    "FW_CLR_PFC05_D", "FW_CLR_PFC1_D", "FW_CLR_PFC2_D"),
                  collapse = "\t"),
           output_file)
close(output_file)

# ------------------------------------------------------------------------------
#  Pull jobs to run, then run 'em!
# ------------------------------------------------------------------------------

wishlist <- read.table(input)
start <- min(start, nrow(wishlist))
end <- min(end, nrow(wishlist))
wishlist <- wishlist[start:end,]

if(any(!(wishlist$baseline %in% c("self", "oracle")))) {
  stop("Invalid reference for some job(s)!")
}
if(any(!(wishlist$partial_info %in% c(0, 1)))) {
  stop("Invalid partial information flag for some job(s)!")
}

output_dir <- file.path("output", "datasets")
results <- NULL

for(i in 1:nrow(wishlist)) {
  cat(paste0("Processing job ", i, " / ", nrow(wishlist), "\n"))

  # Parse data set
  job <- wishlist[i,]
  data <- readRDS(file.path(output_dir, paste0(job$uuid, ".rds")))

  abundances <- data$simulation$abundances
  if(job$partial_info == 1) {
    counts <- data$simulation$observed_counts2
  } else {
    counts <- data$simulation$observed_counts1
  }
  
  med_abs <- mean(rowSums(abundances))
  med_rel <- mean(rowSums(counts))

  n <- nrow(counts)/2

  if(job$type == "relative_abundances") {
    counts_A <- counts[1:n,]
    counts_B <- counts[(n+1):(n*2),]
  } else if(job$type == "scaled_ALDEx2") {
    scaled_counts <- scaled_counts_ALDEx2(counts)
    counts_A <- scaled_counts[1:n,]
    counts_B <- scaled_counts[(n+1):(n*2),]
  } else if(job$type == "scaled_DESeq2") {
    scaled_counts <- scaled_counts_DESeq2(counts, data$simulation$groups)
    counts_A <- scaled_counts[1:n,]
    counts_B <- scaled_counts[(n+1):(n*2),]
  } else if(job$type == "scaled_scran") {
    scaled_counts <- scaled_counts_scran(counts, data$simulation$groups)
    counts_A <- scaled_counts[1:n,]
    counts_B <- scaled_counts[(n+1):(n*2),]
  }

  results_row <- cbind(job, med_abs, med_rel, characterize_dataset(counts_A, counts_B))

  write_delim(results_row, output_fn, delim = "\t", append = TRUE)
}
