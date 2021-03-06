# This is a first pass at a predictive model.

library(ggplot2)
library(gridExtra)
library(betareg)
library(kernlab)

fit_GP <- FALSE

# (1) Parse the simulated data.

label <- "analysis2"
output_file <- file.path("simulated_analyses",label,"results.tsv")

eval_data <- read.table(output_file, header = T, stringsAsFactors = F)
# Note: The latest run includes an insertion into the 6th column of a logical flag that indicates
#       whether or not the NB LRT was used for DE (as opposed to the log-LM + permutations).
#       These labels will have to be updated if using those results.
# eval_data <- eval_data[1:1000,]

# (2) Fit and evaluate a Beta regression model.
# Note: We need to add a tiny value because betareg only allows outcomes (0,1).

addend <- 0.00001
eval_data$fpr <- eval_data$fp / (eval_data$fp + eval_data$tn) + addend
eval_data$fnr <- eval_data$fn / (eval_data$fn + eval_data$tp) + addend

#if(file.exists(file.path("output", "beta_fpr_fit.rds"))) {
#  fit_fpr <- readRDS(file.path("output", "beta_fpr_fit.rds"))
#} else {
  fit_fpr <- betareg(fpr ~ p + prop_da_realized + sf_corr +
                          median_expr_da_quantile + net_dir_da + sparsity + sim_entropy + max_prop_change,
                    data = eval_data, link = "logit")
#}
predicted_fpr <- predict(fit_fpr, eval_data)

#if(file.exists(file.path("output","beta_fnr_fit.rds"))) {
#  fit_fnr <- readRDS(file.path("output", "beta_fnr_fit.rds"))
#} else {
  fit_fnr <- betareg(fnr ~ p + prop_da_realized + sf_corr +
                           median_expr_da_quantile + net_dir_da + sparsity + sim_entropy + max_prop_change,
                     data = eval_data, link = "logit")
#}
predicted_fnr <- predict(fit_fnr, eval_data)

plot_df <- data.frame(predicted_values = predicted_fpr,
                      true_values = eval_data$fpr,
                      p = eval_data$p,
                      proportion_da = eval_data$prop_da_realized,
                      size_factor_correlation = as.factor(eval_data$sf_corr),
                      median_expr_da_quantile = eval_data$median_expr_da_quantile,
                      net_dir_da = eval_data$net_dir_da,
                      sparsity = eval_data$sparsity,
                      sim_entropy = eval_data$sim_entropy,
                      max_prop_change = eval_data$max_prop_change,
                      which = "fpr")
plot_df <- rbind(plot_df,
                 data.frame(predicted_values = predicted_fnr,
                 true_values = eval_data$fnr,
                 p = eval_data$p,
                 proportion_da = eval_data$prop_da_realized,
                 size_factor_correlation = as.factor(eval_data$sf_corr),
                 median_expr_da_quantile = eval_data$median_expr_da_quantile,
                 net_dir_da = eval_data$net_dir_da,
                 sparsity = eval_data$sparsity,
                 sim_entropy = eval_data$sim_entropy,
                 max_prop_change = eval_data$max_prop_change,
                 which = "fnr"))
plot_df$which <- as.factor(plot_df$which)

p <- ggplot(plot_df) +
     geom_point(aes(x = true_values, y = predicted_values, color = p)) +
     facet_wrap(vars(which), nrow = 1, scales = "free")
ggsave(paste0("betareg_model_",label,"_by-p.png"), p, units = "in", dpi = 150, height = 5, width = 10)
p <- ggplot(plot_df) +
     geom_point(aes(x = true_values, y = predicted_values, color = proportion_da)) +
     facet_wrap(vars(which), nrow = 1, scales = "free")
ggsave(paste0("betareg_model_",label,"_by-propDA.png"), p, units = "in", dpi = 150, height = 5, width = 10)

# (3) Fit and evaluate a GP regression model.

if(fit_GP) {
  # x should be N x D
  # y should be N x 1

  # should probably scale these?
  x <- cbind(eval_data$p,
             eval_data$prop_da_realized,
             eval_data$sf_corr,
             eval_data$median_expr_da_quantile,
             eval_data$net_dir_da,
             eval_data$sparsity,
             eval_data$sim_entropy,
             eval_data$max_prop_change)

  y <- eval_data$fpr
  dim(y) <- c(length(y), 1)
  #if(file.exists(file.path("output", "GP_fpr_fit.rds"))) {
  #  fit_fpr <- readRDS(file.path("output", "GP_fpr_fit.rds"))
  #} else {
    fit_fpr <- gausspr(x, y, scaled = TRUE)
  #}
  predicted_fpr <- predict(fit_fpr, x)

  y <- eval_data$fnr
  dim(y) <- c(length(y), 1)
  #if(file.exists(file.path("output", "GP_fnr_fit.rds"))) {
  #  fit_fnr <- readRDS(file.path("output", "GP_fnr_fit.rds"))
  #} else {
    fit_fnr <- gausspr(x, y, scaled = TRUE)
  #}
  # fit can be adjusted via var = ... and kpar = list(sigma = ...)
  predicted_fnr <- predict(fit_fnr, x)

  plot_df <- data.frame(predicted_values = predicted_fpr,
                        true_values = eval_data$fpr,
                        p = eval_data$p,
                        proportion_da = eval_data$prop_da_realized,
                        which = "fpr")
  plot_df <- rbind(plot_df,
                   data.frame(predicted_values = predicted_fnr,
                              true_values = eval_data$fnr,
                              p = eval_data$p,
                              proportion_da = eval_data$prop_da_realized,
                              which = "fnr"))
  plot_df$which <- as.factor(plot_df$which)

  p <- ggplot(plot_df) +
      geom_point(aes(x = true_values, y = predicted_values, color = p)) +
      facet_wrap(vars(which), nrow = 1, scales = "free")
  ggsave(paste0("GPreg_model_",label,"_by-p.png"), p, units = "in", dpi = 150, height = 5, width = 10)
  p <- ggplot(plot_df) +
      geom_point(aes(x = true_values, y = predicted_values, color = proportion_da)) +
      facet_wrap(vars(which), nrow = 1, scales = "free")
  ggsave(paste0("GPreg_model_",label,"_by-propDA.png"), p, units = "in", dpi = 150, height = 5, width = 10)
}
