library(ggplot2)

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

# My hunch about what's happening is we've got really different distributions (i.e. lots of very low
# abundance stuff vs. more moderate abundance stuff) and/or the distribution of genes simulated as
# differentially expressed is very different between groups, causing more "movement" in one case vs.
# the other.

plot_diagnostic <- function(label, filename, components) {
  plot_df <- data.frame(values = c(), type = c())
  for(cc in 1:length(components)) {
    plot_df <- rbind(plot_df, data.frame(values = c(components[[cc]]), type = cc))
  }
  plot_df$type <- as.factor(plot_df$type)

  p <- ggplot(plot_df) +
    geom_density(aes(x = values, color = type)) +
    ggtitle(label)
  ggsave(paste0(filename, ".png"), p, units = "in", dpi = 150, height = 6, width = 6)
}

# (1) Are the baseline distributions different? (NO)
log_abundance_1 <- log(sim1$abundance[sim1$groups == 0,] + 0.5)
log_abundance_2 <- log(sim2$abundance[sim2$groups == 0,] + 0.5)
log_abundance_3 <- log(sim3$abundance[sim3$groups == 0,] + 0.5)
log_abundance_4 <- log(sim4$abundance[sim4$groups == 0,] + 0.5)

plot_diagnostic("baseline log distributions",
                "diagnostic_1",
                list(c(log_abundance_1),
                     c(log_abundance_2),
                     c(log_abundance_3),
                     c(log_abundance_4)))

# (2) Are distributions of (randomly selected) differential expressed genes different at baseline? (NO)

log_abundance_1 <- log(sim1$abundance[sim1$groups == 0,sim1$da_genes] + 0.5)
log_abundance_2 <- log(sim2$abundance[sim2$groups == 0,sim2$da_genes] + 0.5)
log_abundance_3 <- log(sim3$abundance[sim3$groups == 0,sim3$da_genes] + 0.5)
log_abundance_4 <- log(sim4$abundance[sim4$groups == 0,sim4$da_genes] + 0.5)

plot_diagnostic("baseline log distributions (DE genes only)",
                "diagnostic_2",
                list(c(log_abundance_1),
                     c(log_abundance_2),
                     c(log_abundance_3),
                     c(log_abundance_4)))

# (3) Are the distributions of differentially expressed genes different in "treatment" condition? (NO)

log_abundance_1 <- log(sim1$abundance[sim1$groups == 1,sim1$da_genes] + 0.5)
log_abundance_2 <- log(sim2$abundance[sim2$groups == 1,sim2$da_genes] + 0.5)
log_abundance_3 <- log(sim3$abundance[sim3$groups == 1,sim3$da_genes] + 0.5)
log_abundance_4 <- log(sim4$abundance[sim4$groups == 1,sim4$da_genes] + 0.5)

plot_diagnostic("treatment log distributions (DE genes only)",
                "diagnostic_3",
                list(c(log_abundance_1),
                     c(log_abundance_2),
                     c(log_abundance_3),
                     c(log_abundance_4)))


# (4) Is the direction of differential expression different? (MAYBE?)

log_diff_1 <- log(sim1$abundance[sim1$groups == 1,sim1$da_genes] + 0.5) - log(sim1$abundance[sim1$groups == 0,sim1$da_genes] + 0.5)
log_diff_2 <- log(sim2$abundance[sim2$groups == 1,sim2$da_genes] + 0.5) - log(sim2$abundance[sim2$groups == 0,sim2$da_genes] + 0.5)
log_diff_3 <- log(sim3$abundance[sim3$groups == 1,sim3$da_genes] + 0.5) - log(sim3$abundance[sim3$groups == 0,sim3$da_genes] + 0.5)
log_diff_4 <- log(sim4$abundance[sim4$groups == 1,sim4$da_genes] + 0.5) - log(sim4$abundance[sim4$groups == 0,sim4$da_genes] + 0.5)

plot_diagnostic("log differential (DE genes only)",
                "diagnostic_4",
                list(c(log_diff_1),
                     c(log_diff_2),
                     c(log_diff_3),
                     c(log_diff_4)))


