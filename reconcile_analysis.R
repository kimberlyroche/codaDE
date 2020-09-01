# This script just updates the "ledger" that records which files have been processed.

args = commandArgs(trailingOnly=TRUE)
analysis_label <- args[1]

processed_filenames <- list.files(path = file.path("simulated_analyses", analysis_label), pattern = "*.rds")
ledger_file <- file.path("simulated_analyses", analysis_label, "ledger.tsv")
ledger <- data.frame(filename = processed_filenames)
write.table(ledger, file = ledger_file, quote = FALSE, sep = '\t', row.names = FALSE)
