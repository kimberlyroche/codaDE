# PT 1
# pull simulations with high and low false positive rates
results <- read.table("simulated_analyses/analysis2/results.tsv", header = T, stringsAsFactors = FALSE)
results$fpr <- results$fp / (results$fp + results$tn + 0.0001)

# attepmted and realized differential abundance are about 90% correlated

# subset for sanity
p <- 10000
results <- results[results$p == p &
                   results$prop_da_attemped >= 0.8 &
                   results$sf_corr == 0.1,]
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
plot(density(diff1), col = "black", lty = 1, xlim = c(-20000, 20000))
lines(density(diff2), col = "black", lty = 2, type = "l")
lines(density(diff3), col = "red", lty = 1, type = "l")
lines(density(diff4), col = "red", lty = 2, type = "l")
dev.off()

# high error simulations have big net increases in total counts
# are there differences in the true library sizes in the control case? -- effectively no
diff1 <- rowSums(low_error_data[[1]]$data$abundances[1:250,])
diff2 <- rowSums(low_error_data[[2]]$data$abundances[1:250,])
diff3 <- rowSums(high_error_data[[1]]$data$abundances[1:250,])
diff4 <- rowSums(high_error_data[[2]]$data$abundances[1:250,])
png("diagnostic.png")
plot(density(diff1), col = "black", lty = 1, xlim = c(10000, 50000))
lines(density(diff2), col = "black", lty = 2, type = "l")
lines(density(diff3), col = "red", lty = 1, type = "l")
lines(density(diff4), col = "red", lty = 2, type = "l")
dev.off()

# are there differences in the true library sizes in the DIFFERENTIAL case? -- ABSOLUTELY
# the differential expression is NET zero in the low error cases
# NET 1.5x in the high error cases
diff1 <- rowSums(low_error_data[[1]]$data$abundances[251:500,])
diff2 <- rowSums(low_error_data[[2]]$data$abundances[251:500,])
diff3 <- rowSums(high_error_data[[1]]$data$abundances[251:500,])
diff4 <- rowSums(high_error_data[[2]]$data$abundances[251:500,])
png("diagnostic.png")
plot(density(diff1), col = "black", lty = 1, xlim = c(10000, 50000))
lines(density(diff2), col = "black", lty = 2, type = "l")
lines(density(diff3), col = "red", lty = 1, type = "l")
lines(density(diff4), col = "red", lty = 2, type = "l")
dev.off()

# are the selected genes different in their baseline abundances? -- NO
diff1 <- rowSums(low_error_data[[1]]$data$abundances[1:250,low_error_data[[1]]$data$da_genes])
diff2 <- rowSums(low_error_data[[2]]$data$abundances[1:250,low_error_data[[2]]$data$da_genes])
diff3 <- rowSums(high_error_data[[1]]$data$abundances[1:250,high_error_data[[1]]$data$da_genes])
diff4 <- rowSums(high_error_data[[2]]$data$abundances[1:250,high_error_data[[2]]$data$da_genes])
png("diagnostic.png")
plot(density(diff1), col = "black", lty = 1, xlim = c(10000, 50000))
lines(density(diff2), col = "black", lty = 2, type = "l")
lines(density(diff3), col = "red", lty = 1, type = "l")
lines(density(diff4), col = "red", lty = 2, type = "l")
dev.off()

# does the distribution of per-gene changes in mean abundance look different? -- NOT REALLY?
fold1 <- colMeans(low_error_data[[1]]$data$abundances[251:500,low_error_data[[1]]$data$da_genes]) /
            colMeans(low_error_data[[1]]$data$abundances[1:250,low_error_data[[1]]$data$da_genes])
fold2 <- colMeans(low_error_data[[2]]$data$abundances[251:500,low_error_data[[2]]$data$da_genes]) / 
            colMeans(low_error_data[[2]]$data$abundances[1:250,low_error_data[[2]]$data$da_genes])
fold3 <- colMeans(high_error_data[[1]]$data$abundances[251:500,high_error_data[[1]]$data$da_genes]) / 
            colMeans(high_error_data[[1]]$data$abundances[1:250,high_error_data[[1]]$data$da_genes])
fold4 <- colMeans(high_error_data[[2]]$data$abundances[251:500,high_error_data[[2]]$data$da_genes]) / 
            colMeans(high_error_data[[2]]$data$abundances[1:250,high_error_data[[2]]$data$da_genes])
breaks <- seq(2, 100, length.out = 5)
breaks <- c(0, 1, breaks, 1000, 5000, 10000, Inf)
table(cut(fold1, breaks = breaks))
table(cut(fold2, breaks = breaks))
table(cut(fold3, breaks = breaks))
table(cut(fold4, breaks = breaks))

# is this change in total abundances driven by a few genes?
# take only genes below the Kth percentile of differential genes, this should really reduce the difference
qq <- 0.5
diff1 <- rowSums(low_error_data[[1]]$data$abundances[251:500,low_error_data[[1]]$data$da_genes[fold1 < quantile(fold1, probs = c(qq))]])
diff2 <- rowSums(low_error_data[[2]]$data$abundances[251:500,low_error_data[[2]]$data$da_genes[fold2 < quantile(fold2, probs = c(qq))]])
diff3 <- rowSums(high_error_data[[1]]$data$abundances[251:500,high_error_data[[1]]$data$da_genes[fold3 < quantile(fold3, probs = c(qq))]])
diff4 <- rowSums(high_error_data[[2]]$data$abundances[251:500,high_error_data[[2]]$data$da_genes[fold4 < quantile(fold4, probs = c(qq))]])
png("diagnostic_50.png")
plot(density(diff1), col = "black", lty = 1, xlim = c(5, 30000))
lines(density(diff2), col = "black", lty = 2, type = "l")
lines(density(diff3), col = "red", lty = 1, type = "l")
lines(density(diff4), col = "red", lty = 2, type = "l")
dev.off()

qq <- 0.9
diff1 <- rowSums(low_error_data[[1]]$data$abundances[251:500,low_error_data[[1]]$data$da_genes[fold1 < quantile(fold1, probs = c(qq))]])
diff2 <- rowSums(low_error_data[[2]]$data$abundances[251:500,low_error_data[[2]]$data$da_genes[fold2 < quantile(fold2, probs = c(qq))]])
diff3 <- rowSums(high_error_data[[1]]$data$abundances[251:500,high_error_data[[1]]$data$da_genes[fold3 < quantile(fold3, probs = c(qq))]])
diff4 <- rowSums(high_error_data[[2]]$data$abundances[251:500,high_error_data[[2]]$data$da_genes[fold4 < quantile(fold4, probs = c(qq))]])
png("diagnostic_98.png")
plot(density(diff1), col = "black", lty = 1, xlim = c(5, 30000))
lines(density(diff2), col = "black", lty = 2, type = "l")
lines(density(diff3), col = "red", lty = 1, type = "l")
lines(density(diff4), col = "red", lty = 2, type = "l")
dev.off()

qq <- 1
diff1 <- rowSums(low_error_data[[1]]$data$abundances[251:500,low_error_data[[1]]$data$da_genes[fold1 < quantile(fold1, probs = c(qq))]])
diff2 <- rowSums(low_error_data[[2]]$data$abundances[251:500,low_error_data[[2]]$data$da_genes[fold2 < quantile(fold2, probs = c(qq))]])
diff3 <- rowSums(high_error_data[[1]]$data$abundances[251:500,high_error_data[[1]]$data$da_genes[fold3 < quantile(fold3, probs = c(qq))]])
diff4 <- rowSums(high_error_data[[2]]$data$abundances[251:500,high_error_data[[2]]$data$da_genes[fold4 < quantile(fold4, probs = c(qq))]])
png("diagnostic_100.png")
plot(density(diff1), col = "black", lty = 1, xlim = c(5, 30000))
lines(density(diff2), col = "black", lty = 2, type = "l")
lines(density(diff3), col = "red", lty = 1, type = "l")
lines(density(diff4), col = "red", lty = 2, type = "l")
dev.off()








# PT 2
# Visualize some erroneous DA calls
# data <- readRDS(paste0("simulated_data/",high_error_fns[1]))

# non_da_genes <- !(1:ncol(data$data$abundances) %in% data$data$da_genes)
# evaluate_features <- apply(data$data$observed_counts, 2, function(x) mean(x) > 3)
# plot_genes <- which(non_da_genes & evaluate_features)

# non_da_gene <- sample(plot_genes)[1]

# png("before.png")
# plot(data$data$abundances[,non_da_gene])
# dev.off()

# png("after.png")
# plot(data$data$observed_counts[,non_da_gene])
# dev.off()

