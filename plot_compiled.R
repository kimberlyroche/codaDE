library(tidyverse)
library(codaDE)
library(ggROC)
library(pROC)

setwd("C:/Users/kimbe/Documents/codaDE")

p <- 100

data <- readRDS(file.path("output", paste0("simresults_p",p,"_simulated_all.rds")))

# Peek at realized fold changes
# hist(data$delta_mean_v2)

# ------------------------------------------------------------------------------
#   ROC curve(ish) version
# ------------------------------------------------------------------------------

plot_data <- data %>%
  select(delta_mean_v2, rate, rate_type, method) %>%
  filter(method != "spike_in") %>%
  pivot_wider(names_from = rate_type, values_from = rate)
plot_data$method <- factor(plot_data$method)
levels(plot_data$method) <- c("NB GLM", "DESeq2", "MAST", "scran")

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
ggsave(file.path("output", "images", save_filename <- paste0("DE_p",p,"_all_models.png")),
       units = "in",
       height = 5,
       width = 5)

# ------------------------------------------------------------------------------
#   Partial abundance version (with NB GLM only)
# ------------------------------------------------------------------------------

partial_plot <- TRUE

if(partial_plot) {
  plot_data <- data %>%
    select(delta_mean_v2, rate, rate_type, method) %>%
    filter(method %in% c("baseline", "spike_in")) %>%
    pivot_wider(names_from = rate_type, values_from = rate)
  plot_data$method <- factor(plot_data$method)
  levels(plot_data$method) <- c("NB GLM", "NB GLM (partial totals)")

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
  ggsave(file.path("output", "images", paste0("DE_p",p,"_all_models_partial.png")),
         units = "in",
         height = 2.5,
         width = 5)
}

# ------------------------------------------------------------------------------
#   DESeq2 performance (FPR) as a function of log fold difference in means
# ------------------------------------------------------------------------------

lfc_plot <- TRUE

if(lfc_plot) {
  data <- cbind(p = 100, readRDS(file.path("output", "simresults_p100_simulated_all.rds")))
  for(p in c(1000, 5000, 15000)) {
    data <- rbind(data, cbind(p = p, readRDS(file.path("output", paste0("simresults_p",p,"_simulated_all.rds")))))
  }
  data <- data %>%
    filter(method == "DESeq2") %>%
    filter(p > 100) %>%
    filter(rate_type == "fpr") %>%
    filter(delta_mean_v2 <= 10) %>%
    select(p, delta_mean_v2, rate, method)
  data$p <- factor(data$p)
  
  ggplot(data, aes(x = delta_mean_v2, y = rate, color = p)) + # plot an ROC curve
    geom_point(size = 2, alpha = 0.5) +
    geom_smooth(method = "loess") +
    scale_color_manual(values = generate_highcontrast_palette(length(levels(data$p)))) #+
    # facet_wrap(. ~ p)
}

# ------------------------------------------------------------------------------
#   DESeq2 performance (FPR) as a function of log fold difference in means
# ------------------------------------------------------------------------------

data <- plot_data %>%
  filter(method %in% c("baseline", "spike_in")) %>%
  filter(rate_type == "fpr") %>%
  filter(cor_totals >= 0) %>%
  select(cor_totals, rate, method)

ggplot(data, aes(x = cor_totals, y = rate)) + # plot an ROC curve
  geom_point(size = 2, alpha = 0.5)

# ------------------------------------------------------------------------------
#   FPR/TPR version
# ------------------------------------------------------------------------------

# rt <- "fpr"
# sub_data <- data %>%
#   filter(rate_type == rt) %>%
#   filter(delta_mean_v2 < 10) %>%
#   filter(method %in% c("baseline", "DESeq2", "MAST", "scran"))
# 
# sub_data$method <- factor(sub_data$method)
# levels(sub_data$method) = c("NB GLM",
#                             "DESeq2",
#                             "MAST",
#                             "scran")
# 
# pl <- ggplot(sub_data, aes(x = log(delta_mean_v2), y = rate, color = method)) +
#   scale_color_manual(values = cpalette) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(alpha = 0.2) +
#   labs(color = "method")
# if(rt == "tpr") {
#   pl <- pl +
#     ylim(c(0, 1))
# } else {
#   pl <- pl +
#     ylim(c(0, 0.45))
# }
# pl <- pl +
#   xlab("log fold change in abundance between conditions")
# if(tpr) {
#   pl <- pl +
#     ylab("true positive rate")
# } else {
#   pl <- pl +
#     ylab("false positive rate")
# }
# 
# show(pl)
# 
# ggsave(file.path("output", "images", paste0("DE_p",p,"_all_models_",rt,".png")),
#        pl,
#        dpi = 500,
#        units = "in",
#        # height = 6.5,
#        # width = 9)
#        height = 4,
#        width = 6)



