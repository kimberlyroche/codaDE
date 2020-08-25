# This is a first pass at a predictive model.

library(ggplot2)
library(gridExtra)
library(betareg)
library(kernlab)

# (1) Parse the simulated data.

output_file <- file.path("results.txt")

data <- read.table(output_file, header = T)
# Note: The latest run includes an insertion into the 6th column of a logical flag that indicates
#       whether or not the NB LRT was used for DE (as opposed to the log-LM + permutations).
#       These labels will have to be updated if using those results.

# (2) Fit and evaluate a Beta regression model.
# Note: We need to add a tiny value because betareg only allows outcomes (0,1).

eval_data <- data
eval_data <- eval_data[eval_data$prop_de > 0,]
eval_data <- eval_data[eval_data$prop_de < 1,]
eval_data <- eval_data[eval_data$p >= 1000,]

addend <- 0.00001
eval_data$fpr <- eval_data$fp / (eval_data$fp + eval_data$tn) + addend
eval_data$fnr <- eval_data$fn / (eval_data$fn + eval_data$tp) + addend

fit_fpr <- betareg(fpr ~ p + prop_de + sf_corr + expr_de + dir_de + sparsity + entropy,
                   data = eval_data, link = "log")
predicted_fpr <- predict(fit_fpr, eval_data)

fit_fnr <- betareg(fnr ~ p + prop_de + sf_corr + expr_de + dir_de + sparsity + entropy,
                   data = eval_data, link = "log")
predicted_fnr <- predict(fit_fnr, eval_data)

plot_df <- data.frame(predicted_values = predicted_fpr,
                      true_values = eval_data$fpr,
                      p = as.factor(eval_data$p),
                      prop_de = as.factor(eval_data$prop_de),
                      sf_corr = as.factor(eval_data$sf_corr),
                      expr_de = eval_data$expr_de,
                      dir_de = eval_data$dir_de,
                      sparsity = eval_data$sparsity,
                      entropy = eval_data$entropy,
                      which = "fpr")
plot_df <- rbind(plot_df,
                 data.frame(predicted_values = predicted_fnr,
                 true_values = eval_data$fnr,
                 p = as.factor(eval_data$p),
                 prop_de = as.factor(eval_data$prop_de),
                 sf_corr = as.factor(eval_data$sf_corr),
                 expr_de = eval_data$expr_de,
                 dir_de = eval_data$dir_de,
                 sparsity = eval_data$sparsity,
                 entropy = eval_data$entropy,
                 which = "fnr"))
plot_df$which <- as.factor(plot_df$which)

p <- ggplot(plot_df) +
     geom_point(aes(x = true_values, y = predicted_values, color = p)) +
     facet_wrap(vars(which), nrow = 1, scales = "free")
ggsave("betareg_model.png", p, units = "in", dpi = 150, height = 5, width = 10)
