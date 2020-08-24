# This is a first pass at a predictive model.

library(ggplot2)
library(gridExtra)
library(betareg)
library(kernlab)

# (1) Parse the simulated data.

output_file <- file.path("output", "results_singlecellRNAseq.tsv")
data <- read.table(output_file, header = F)
colnames(data) <- c("run_label", "p", "prop_de", "sfcorr", "measure", "DE_NB",
                    "min_abundance", "tp", "fp", "tn", "fn", "feature_evaluated")
# Note: The latest run includes an insertion into the 6th column of a logical flag that indicates
#       whether or not the NB LRT was used for DE (as opposed to the log-LM + permutations).
#       These labels will have to be updated if using those results.

# (2) Calculate and append true positive rate, etc.
data$tpr <- data$tp / (data$tp + data$fn)
data$fpr <- data$fp / (data$fp + data$tn)
data$tnr <- data$tn / (data$fp + data$tn)
data$fnr <- data$fn / (data$tp + data$fn)

# (3) Plot the data. Does anything have an obvious linear or log-linear relationship with 
#     the outcome?

eval_data <- data[data$DE_NB == TRUE &
                  data$prop_de > 0,]

# outcome isn't really linear in any of the predictors
p1 <- ggplot(eval_data) +
       geom_point(aes(x = p, y = fpr))
p2 <- ggplot(eval_data) +
       geom_point(aes(x = p, y = fnr))
p3 <- ggplot(eval_data) +
       geom_point(aes(x = prop_de, y = fpr))
p4 <- ggplot(eval_data) +
       geom_point(aes(x = prop_de, y = fnr))
p5 <- ggplot(eval_data) +
       geom_point(aes(x = sfcorr, y = fpr))
p6 <- ggplot(eval_data) +
       geom_point(aes(x = sfcorr, y = fnr))
p <- grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)
ggsave("predictor_correlation.png", p)

# (4) Fit and evaluate a Beta regression model.
# Note: We need to add a tiny value because betareg only allows outcomes (0,1).

addend <- 0.00001
eval_data$fpr <- eval_data$fpr + addend
eval_data$fnr <- eval_data$fnr + addend

fit_fpr <- betareg(fpr ~ p + prop_de + sfcorr, data = eval_data, link = "log")
predicted_fpr <- predict(fit_fpr, eval_data)

fit_fnr <- betareg(fnr ~ p + prop_de + sfcorr, data = eval_data, link = "log")
predicted_fnr <- predict(fit_fnr, eval_data)

plot_df <- data.frame(predicted_values = predicted_fpr, true_values = eval_data$fpr, which = "fpr")
plot_df <- rbind(plot_df,
                 data.frame(predicted_values = predicted_fnr, true_values = eval_data$fnr, which = "fnr"))
plot_df$which <- as.factor(plot_df$which)

p <- ggplot(plot_df) +
     geom_point(aes(x = true_values, y = predicted_values)) +
     facet_wrap(vars(which), nrow = 1, scales = "free")
ggsave("betareg_model.png", p, units = "in", dpi = 150, height = 5, width = 10)

# Hopefully adding gene number improves this a lot!

# (5) Fit and evaluate a GP regression model.
# Note: Fit can be adjusted via var = ... and kpar = list(sigma = ...).

# x should be N x D
# y should be N x 1

x <- cbind(eval_data$p, eval_data$prop_de, eval_data$sfcorr)
y <- eval_data$fpr

dim(y) <- c(length(y), 1)
fit_fpr <- gausspr(x, y, scaled = TRUE)
predicted_fpr <- predict(fit_fpr, x)

y <- eval_data$fnr
dim(y) <- c(length(y), 1)
fit_fnr <- gausspr(x, y, scaled = TRUE)
# fit can be adjusted via var = ... and kpar = list(sigma = ...)
predicted_fnr <- predict(fit_fnr, x)

plot_df <- data.frame(predicted_values = predicted_fpr, true_values = eval_data$fpr, which = "fpr")
plot_df <- rbind(plot_df,
                 data.frame(predicted_values = predicted_fnr, true_values = eval_data$fnr, which = "fnr"))
plot_df$which <- as.factor(plot_df$which)

p <- ggplot(plot_df) +
     geom_point(aes(x = true_values, y = predicted_values)) +
     facet_wrap(vars(which), nrow = 1, scales = "free")
ggsave("GPreg_model.png", p, units = "in", dpi = 150, height = 5, width = 10)



