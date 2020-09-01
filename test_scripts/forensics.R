# PT 1
# pull simulations with high and low false positive rates
results <- read.table("simulated_data/results_filter0.tsv", header = T, stringsAsFactors = FALSE)
results$fpr <- results$fp / (results$fp + results$tn)

# subset for sanity
results <- results[results$p == 10000 &
                   results$proportion_da == 0.8 &
                   results$size_factor_correlation == 0.1,]
head(results)
round(min(results$fpr), 3)
round(max(results$fpr), 3)
results <- results[order(results$fpr),]

# now we have simulations that were created under the same conditions but have very different
# characteristic false positive rates!
low_error_fns <- results[1:2,]$filename
high_error_fns <- results[(nrow(results) - 1):nrow(results),]$filename

low_error_data <- list(readRDS(file.path("simulated_data", low_error_fns[1])),
                       readRDS(file.path("simulated_data", low_error_fns[2])))

high_error_data <- list(readRDS(file.path("simulated_data", high_error_fns[1])),
                        readRDS(file.path("simulated_data", high_error_fns[2])))

# lots of variation in error is associated with LOW/NO information about original total abundances
# is there a systematic difference in the resampling?
# (a) compare total abundances (before and after resampling) between the low and high error cases
diff1 <- rowSums(low_error_data[[1]]$data$abundances) - rowSums(low_error_data[[1]]$data$observed_counts)
diff2 <- rowSums(low_error_data[[2]]$data$abundances) - rowSums(low_error_data[[2]]$data$observed_counts)
diff3 <- rowSums(high_error_data[[1]]$data$abundances) - rowSums(high_error_data[[1]]$data$observed_counts)
diff4 <- rowSums(high_error_data[[2]]$data$abundances) - rowSums(high_error_data[[2]]$data$observed_counts)
png("diagnostic1.png")
plot(density(diff1), col = "black", lty = 1)
lines(density(diff2), col = "black", lty = 2, type = "l")
lines(density(diff3), col = "red", lty = 1, type = "l")
lines(density(diff4), col = "red", lty = 2, type = "l")
dev.off()


# PT 2
# Visualize some erroneous DA calls
data <- readRDS("simulated_data/0170698d-1aea-4142-b320-08971893a0ae.rds")

non_da_genes <- !(1:ncol(data$data$abundances) %in% data$data$da_genes)
evaluate_features <- apply(data$data$observed_counts, 2, function(x) mean(x) > 3)
plot_genes <- which(non_da_genes & evaluate_features)

non_da_gene <- sample(plot_genes)[1]

png("before.png")
plot(data$data$abundances[,non_da_gene])
dev.off()

png("after.png")
plot(data$data$observed_counts[,non_da_gene])
dev.off()

