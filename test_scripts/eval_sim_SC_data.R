simdata <- list()
for(i in 1:20) {
  simdata[[i]] <- readRDS(paste0("sim_dataset_",i,".rds"))
}

# show error
tpr_vector <- c()
fpr_vector <- c()
for(i in 1:length(simdata)) {
  dataset <- simdata[[i]]
  tpr <- dataset$TP / (dataset$TP + dataset$FN)
  fpr <- dataset$FP / (dataset$FP + dataset$TN)
  tpr_vector <- c(tpr_vector, tpr)
  fpr_vector <- c(fpr_vector, fpr)
}

# what's different about simulations that have very different TPR's?
# sim1 <- simdata[[order(tpr_vector)[20]]]$data # "good"
# sim2 <- simdata[[order(tpr_vector)[19]]]$data # "good"
# sim3 <- simdata[[order(tpr_vector)[2]]]$data # "bad"
# sim4 <- simdata[[order(tpr_vector)[1]]]$data # "bad"

# ... or FPR's?
sim1 <- simdata[[order(fpr_vector)[1]]]$data # "good"
sim2 <- simdata[[order(fpr_vector)[2]]]$data # "good"
sim3 <- simdata[[order(fpr_vector)[19]]]$data # "bad"
sim4 <- simdata[[order(fpr_vector)[20]]]$data # "bad"

# (1) do the differentially expressed genes have a different mean log abundance?

da_mean_log_ab1 <- colMeans(log(sim1$observed_counts[1:250,sim1$da_genes] + 0.5))
da_mean_log_ab2 <- colMeans(log(sim2$observed_counts[1:250,sim2$da_genes] + 0.5))
da_mean_log_ab3 <- colMeans(log(sim3$observed_counts[1:250,sim3$da_genes] + 0.5))
da_mean_log_ab4 <- colMeans(log(sim4$observed_counts[1:250,sim4$da_genes] + 0.5))
png("diagnostic_1a.png")
plot(density(da_mean_log_ab1))
lines(density(da_mean_log_ab2), col = "black", lty = 2)
lines(density(da_mean_log_ab3), col = "red")
lines(density(da_mean_log_ab4), col = "red", lty = 2)
dev.off()

n_genes <- ncol(sim1$abundances)

nonda_mean_log_ab1 <- colMeans(log(sim1$observed_counts[1:250,setdiff(1:n_genes, sim1$da_genes)] + 0.5))
nonda_mean_log_ab2 <- colMeans(log(sim2$observed_counts[1:250,setdiff(1:n_genes, sim2$da_genes)] + 0.5))
nonda_mean_log_ab3 <- colMeans(log(sim3$observed_counts[1:250,setdiff(1:n_genes, sim3$da_genes)] + 0.5))
nonda_mean_log_ab4 <- colMeans(log(sim4$observed_counts[1:250,setdiff(1:n_genes, sim4$da_genes)] + 0.5))
png("diagnostic_1b.png")
plot(density(nonda_mean_log_ab1))
lines(density(nonda_mean_log_ab2), col = "black", lty = 2)
lines(density(nonda_mean_log_ab3), col = "red")
lines(density(nonda_mean_log_ab4), col = "red", lty = 2)
dev.off()

# (2) is the baseline distribution different (i.e. in terms of uniformity)?

mean_log_ab1 <- colMeans(log(sim1$observed_counts[1:250,] + 0.5))
mean_log_ab2 <- colMeans(log(sim2$observed_counts[1:250,] + 0.5))
mean_log_ab3 <- colMeans(log(sim3$observed_counts[1:250,] + 0.5))
mean_log_ab4 <- colMeans(log(sim4$observed_counts[1:250,] + 0.5))
png("diagnostic_2.png")
plot(density(mean_log_ab1))
lines(density(mean_log_ab2), col = "black", lty = 2)
lines(density(mean_log_ab3), col = "red")
lines(density(mean_log_ab4), col = "red", lty = 2)
dev.off()

# (3) is the direction of differential expression skewed positive or negative?

log_da_counts1a <- colMeans(log(sim1$observed_counts[1:250,sim1$da_genes] + 0.5))
log_da_counts1b <- colMeans(log(sim1$observed_counts[251:500,sim1$da_genes] + 0.5))
direction1 <- log_da_counts1a - log_da_counts1b

log_da_counts2a <- colMeans(log(sim2$observed_counts[1:250,sim2$da_genes] + 0.5))
log_da_counts2b <- colMeans(log(sim2$observed_counts[251:500,sim2$da_genes] + 0.5))
direction2 <- log_da_counts2a - log_da_counts2b

log_da_counts3a <- colMeans(log(sim3$observed_counts[1:250,sim3$da_genes] + 0.5))
log_da_counts3b <- colMeans(log(sim3$observed_counts[251:500,sim3$da_genes] + 0.5))
direction3 <- log_da_counts3a - log_da_counts3b

log_da_counts4a <- colMeans(log(sim4$observed_counts[1:250,sim4$da_genes] + 0.5))
log_da_counts4b <- colMeans(log(sim4$observed_counts[251:500,sim4$da_genes] + 0.5))
direction4 <- log_da_counts4a - log_da_counts4b

png("diagnostic_3.png")
plot(density(direction1))
lines(density(direction2), col = "black", lty = 2)
lines(density(direction3), col = "red")
lines(density(direction4), col = "red", lty = 2)
dev.off()

# (4) is there a net difference in total abundance?

total_ab1a <- rowSums(sim1$abundances[1:250,])
total_ab1b <- rowSums(sim1$abundances[251:500,])
difference1 <- total_ab1b - total_ab1a

total_ab2a <- rowSums(sim2$abundances[1:250,])
total_ab2b <- rowSums(sim2$abundances[251:500,])
difference2 <- total_ab2b - total_ab2a

total_ab3a <- rowSums(sim3$abundances[1:250,])
total_ab3b <- rowSums(sim3$abundances[251:500,])
difference3 <- total_ab3b - total_ab3a

total_ab4a <- rowSums(sim4$abundances[1:250,])
total_ab4b <- rowSums(sim4$abundances[251:500,])
difference4 <- total_ab4b - total_ab4a

png("diagnostic_4.png")
plot(density(difference1))
lines(density(difference2), col = "black", lty = 2)
lines(density(difference3), col = "red")
lines(density(difference4), col = "red", lty = 2)
dev.off()



