library(ggplot2)

for(cond in c(1,2)) {
  pattern_str <- paste0("results_", cond, "_.*?\\.rds")
  results_dir <- "output"
  results_files <- list.files(path = results_dir, pattern = pattern_str, full.names = TRUE, recursive = FALSE)
  full_results <- NULL
  for(rf in results_files) {
    if(is.null(full_results)) {
      full_results <- readRDS(rf)
    } else {
      full_results <- rbind(full_results, readRDS(rf))
    }
  }

  if(cond == 1) {
    p <- ggplot(full_results, aes(x = as.factor(prop_de), y = FPR)) +
      geom_boxplot() +
      ylim(0, 0.3) +
      xlab("approx. proportion DE genes")
    ggsave(paste0("FPR_", cond, ".png"), p, units = "in", dpi = 100, height = 5, width = 5)
    p <- ggplot(full_results, aes(x = as.factor(prop_de), y = TPR)) +
      geom_boxplot() +
      ylim(0.75, 1) +
      xlab("approx. proportion DE genes")
    ggsave(paste0("TPR_", cond, ".png"), p, units = "in", dpi = 100, height = 5, width = 5)
  } else {
    p <- ggplot(full_results, aes(x = as.factor(sf_corr), y = FPR)) +
      geom_boxplot() +
      ylim(0, 0.3) +
      xlab("approx. correlation true/observed library sizes")
    ggsave(paste0("FPR_", cond, ".png"), p, units = "in", dpi = 100, height = 5, width = 5)
    p <- ggplot(full_results, aes(x = as.factor(sf_corr), y = TPR)) +
      geom_boxplot() +
      ylim(0.75, 1) +
      xlab("approx. correlation true/observed library sizes")
    ggsave(paste0("TPR_", cond, ".png"), p, units = "in", dpi = 100, height = 5, width = 5)
  }
}


