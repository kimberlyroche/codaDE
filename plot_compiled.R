library(tidyverse)
library(codaDE)
library(ggROC)
library(pROC)

setwd("C:/Users/kimbe/Documents/codaDE")

p <- 5000

data <- readRDS(file.path("output", paste0("simresults_p",p,"_simulated_all.rds")))

# Peek at realized fold changes
# hist(data$delta_mean_v2)

# ------------------------------------------------------------------------------
#   ROC curve version
# ------------------------------------------------------------------------------

spikeins <- FALSE

if(spikeins) {
  plot_data <- data %>%
    select(delta_mean_v2, rate, rate_type, method) %>%
    filter(method %in% c("baseline", "spike_in")) %>%
    pivot_wider(names_from = rate_type, values_from = rate)
  plot_data$method <- factor(plot_data$method)
  levels(plot_data$method) <- c("NB GLM", "NB GLM (partial totals)")
} else {
  plot_data <- data %>%
    select(delta_mean_v2, rate, rate_type, method) %>%
    filter(method != "spike_in") %>%
    pivot_wider(names_from = rate_type, values_from = rate)
  plot_data$method <- factor(plot_data$method)
  levels(plot_data$method) <- c("NB GLM", "DESeq2", "MAST", "scran")
}

head(plot_data)

if(!exists("cpalette")) {
  cpalette <- generate_highcontrast_palette(100)
  cpalette <- sample(cpalette, size = length(levels(plot_data$method)))
}

ggplot(plot_data, aes(x = fpr, y = tpr, color = method)) + # plot an ROC curve
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = cpalette) +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  xlab("FPR") +
  ylab("TPR") +
  facet_wrap(. ~ method) +
  theme(legend.position = "none")

if(spikeins) {
  ggsave(file.path("output", "images", paste0("DE_p",p,"_all_models_partial.png")),
         units = "in",
         height = 2.5,
         width = 5)
} else {
  ggsave(file.path("output", "images", save_filename <- paste0("DE_p",p,"_all_models.png")),
         units = "in",
         height = 5,
         width = 5)
}

# ------------------------------------------------------------------------------
#   FPR/TPR version
# ------------------------------------------------------------------------------

if(FALSE) {
  rt <- "fpr"
  sub_data <- data %>%
    filter(rate_type == rt) %>%
    filter(delta_mean_v2 < 10) %>%
    filter(method %in% c("baseline", "DESeq2", "MAST", "scran"))
  
  sub_data$method <- factor(sub_data$method)
  levels(sub_data$method) = c("NB GLM",
                              "DESeq2",
                              "MAST",
                              "scran")
  
  pl <- ggplot(sub_data, aes(x = log(delta_mean_v2), y = rate, color = method)) +
    scale_color_manual(values = cpalette) +
    geom_point(alpha = 0.5) +
    geom_smooth(alpha = 0.2) +
    labs(color = "method")
  if(rt == "tpr") {
    pl <- pl +
      ylim(c(0, 1))
  } else {
    pl <- pl +
      ylim(c(0, 0.45))
  }
  pl <- pl +
    xlab("log fold change in abundance between conditions")
  if(tpr) {
    pl <- pl +
      ylab("true positive rate")
  } else {
    pl <- pl +
      ylab("false positive rate")
  }
  
  show(pl)
  
  ggsave(file.path("output", "images", paste0("DE_p",p,"_all_models_",rt,".png")),
         pl,
         dpi = 500,
         units = "in",
         # height = 6.5,
         # width = 9)
         height = 4,
         width = 6)
}


