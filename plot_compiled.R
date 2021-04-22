library(tidyverse)
library(codaDE)
library(ggROC)
library(pROC)

setwd("C:/Users/kimbe/Documents/codaDE")

# Select global palette for 5 methods
# palette <- generate_highcontrast_palette(1000)
# palette <- sample(palette, size = 5)
# points <- data.frame(x = rnorm(5), y = rnorm(5), label = factor(1:5))
# ggplot(points, aes(x = x, y = y, color = label)) +
#   geom_point(size = 10) +
#   scale_color_manual(values = palette)

# Previous
# palette <- c("#46A06B", "#B95D6E", "#EF82BB", "#755A7F", "#E3C012")

palette <- c("#46A06B", "#FF5733", "#EF82BB", "#755A7F", "#E3C012", "#B95D6E")

# ------------------------------------------------------------------------------
#   ROC curve(ish) version
# ------------------------------------------------------------------------------

p <- 1000
data <- readRDS(file.path("output", paste0("simresults_p",p,"_simulated_all.rds")))

plot_data <- data %>%
  select(delta_mean_v2, rate, rate_type, method) %>%
  filter(!(method %in% c("spike_in", "ALDEx2"))) %>%
  pivot_wider(names_from = rate_type, values_from = rate)
plot_data$method <- factor(plot_data$method)
levels(plot_data$method) <- c("NB GLM", "DESeq2", "MAST", "scran")

ggplot(plot_data, aes(x = fpr, y = tpr, color = method)) + # plot an ROC curve
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = palette) +
  # xlim(c(0,1)) +
  # ylim(c(0,1)) +
  xlim(c(0,0.75)) +
  ylim(c(0.25,1)) +
  xlab("FPR") +
  ylab("TPR") +
  facet_wrap(. ~ method) +
  theme(legend.position = "none")
ggsave(file.path("output", "images", save_filename <- paste0("DE_p",p,"_all_models.png")),
       units = "in",
       height = 5,
       width = 5)

# ------------------------------------------------------------------------------
#   Partial abundance version (with NB GLM only) TBD w/ data
# ------------------------------------------------------------------------------

partial_plot <- TRUE

if(partial_plot) {
  p <- 100
  data <- readRDS(file.path("output", paste0("simresults_p",p,"_simulated_all.rds")))
  data <- data %>%
    filter(method %in% c("baseline", "spike_in")) %>%
    filter(rate_type == "fpr") %>%
    filter(delta_mean_v2 > 2)
  
  plot_data <- data.frame(cor_totals = c(), delta_fpr = c())
  for(i in seq(from = 1, to = (nrow(data)-1), by = 2)) {
    delta_fpr <- data$rate[i+1] - data$rate[i]
    cor_totals <- data$cor_totals[i+1]
    plot_data <- rbind(plot_data,
                       data.frame(cor_totals = cor_totals, delta_fpr = delta_fpr))
  }
  ggplot(plot_data, aes(x = cor_totals, y = delta_fpr)) +
    geom_point(size = 1) +
    geom_smooth(method = "lm", color = "#888888", alpha = 0.3)
  ggsave(file.path("output", "images", "partial_total_performance.png"),
         units = "in",
         height = 3,
         width = 3)
  
  # ROC-like plots for NB GLM and NB GLM on partially faithful totals
  data <- readRDS(file.path("output", paste0("simresults_p",p,"_simulated_all.rds")))

  plot_data <- data %>%
    select(delta_mean_v2, rate, rate_type, method, cor_totals) %>%
    # filter(method %in% c("baseline", "spike_in", "ALDEx2")) %>%
    filter(method %in% c("ALDEx2", "DESeq2")) %>%
    pivot_wider(names_from = rate_type, values_from = rate)
  plot_data$method <- factor(plot_data$method)
  # levels(plot_data$method) <- c("NB GLM", "NB GLM (partial totals)", "ALDEx2")
  levels(plot_data$method) <- c("ALDEx2", "DESeq2")

  # binary_palette <- palette[c(1,5,6)]
  binary_palette <- palette[c(6,2)]

  ggplot(plot_data, aes(x = fpr, y = tpr, color = method)) +
    geom_point(size = 2, alpha = 0.5) +
    scale_color_manual(values = binary_palette) +
    # xlim(c(0,1)) +
    # ylim(c(0,1)) +
    xlim(c(0,0.75)) +
    ylim(c(0,1)) +
    xlab("FPR") +
    ylab("TPR") +
    facet_wrap(. ~ method) +
    theme(legend.position = "none")
  ggsave(file.path("output",
                   "images",
                   paste0("partial_total_performance_wALDEx2_p",p,".png")),
         units = "in",
         height = 2.5,
         width = 5)
}

# ------------------------------------------------------------------------------
#   DESeq2 performance (FPR) as a function of log fold difference in means
# ------------------------------------------------------------------------------

lfc_plot <- TRUE

if(lfc_plot) {
  unary_palette <- palette[2]
  
  pp <- c(1000, 5000, 15000)
  data <- cbind(p = 1000, readRDS(file.path("output", paste0("simresults_p",pp[1],"_simulated_all.rds"))))
  for(p in pp[2:length(pp)]) {
    data <- rbind(data, cbind(p = p, readRDS(file.path("output", paste0("simresults_p",p,"_simulated_all.rds")))))
  }
  data <- data %>%
    filter(method == "DESeq2") %>%
    filter(p > 100) %>%
    filter(rate_type == "fpr") %>%
    filter(delta_mean_v2 <= 10) %>%
    select(p, delta_mean_v2, rate, rate_type, method)
  data$p <- factor(data$p)
  
  ggplot(data, aes(x = log(delta_mean_v2), y = rate, color = rate_type)) + # plot an ROC curve
    geom_point(size = 1, alpha = 0.25) +
    geom_smooth(method = "loess", alpha = 0.5) +
    scale_color_manual(values = unary_palette) +
    theme(legend.position = "none") +
    ylim(c(0,0.15))
  ggsave(file.path("output", "images", "DESeq2_performance_curves.png"),
         units = "in",
         height = 3,
         width = 3)
}

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



