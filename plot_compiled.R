source("path_fix.R")

library(tidyverse)
library(codaDE)
library(ggROC)
library(pROC)
library(gridExtra)

source("ggplot_fix.R")

setwd("C:/Users/kimbe/Documents/codaDE")

# Select global palette for 5 methods
# palette <- generate_highcontrast_palette(1000)
# palette <- sample(palette, size = 5)
# points <- data.frame(x = rnorm(5), y = rnorm(5), label = factor(1:5))
# ggplot(points, aes(x = x, y = y, color = label)) +
#   geom_point(size = 10) +
#   scale_color_manual(values = palette)

# Previous
# palette <- generate_highcontrast_palette(7)
# palette <- c("#46A06B", "#B95D6E", "#EF82BB", "#755A7F", "#E3C012")

palette <- c("#46A06B", "#FF5733", "#EF82BB", "#755A7F", "#E3C012", "#B95D6E")

eliminate_outliers <- function(df) {
  # Eliminate outlier simulations as above
  # These are ones with very small or very large differences between 
  #   conditions
  log_diff <- log(df$delta_mean_v2)
  mean_log_diff <- mean(log_diff)
  limit <- 3*sd(log_diff)
  idx_retain <- which(log_diff > mean_log_diff - limit & log_diff < mean_log_diff + limit)
  return(list(df = df[idx_retain,],
              log_diff = log_diff,
              mean_log_diff = mean_log_diff))
}

# ------------------------------------------------------------------------------
#   TPR x FPR plots
# ------------------------------------------------------------------------------

# Note: There are edgeR results in here but I think they are wrong. (They
# perfectly resemble the NB GLM results.)

p_all <- c(100, 1000, 5000)
corrp_all <- c(0, 0.5)
lfc_data <- data.frame(lfc = c(), p = c(), corrp = c())
for(p in p_all) {
  for(corrp in corrp_all) {
    data <- readRDS(file.path("output", paste0("simresults_p",p,"_corrp",corrp,"_all.rds")))

    plot_data <- data %>%
      select(delta_mean_v1, delta_mean_v2, rate, rate_type, method) %>%
      filter(method %in% c("baseline", "DESeq2", "MAST", "scran", "ALDEx2")) %>%
      pivot_wider(names_from = rate_type, values_from = rate)
    # Assign factor order...
    plot_data$method <- factor(plot_data$method, levels = c("baseline",
                                                            "MAST",
                                                            "DESeq2",
                                                            "scran",
                                                            "ALDEx2"))
    # ...then rename
    levels(plot_data$method) <- c("NB GLM",
                                  "MAST",
                                  "DESeq2",
                                  "scran",
                                  "ALDEx2")
    
    
    plot_data <- eliminate_outliers(plot_data)
    log_diff <- plot_data$log_diff
    mean_log_diff <- plot_data$mean_log_diff
    plot_data <- plot_data$df
    
    lfc_data <- rbind(lfc_data,
                      data.frame(lfc = log_diff, p = p, corrp = corrp))
    
    ggplot(plot_data, aes(x = fpr, y = tpr, color = method)) +
      geom_point(size = 1, alpha = 0.5) +
      scale_color_manual(values = palette) +
      xlim(c(0,0.75)) +
      ylim(c(0.25,1)) +
      xlab("1 - specificity") +
      ylab("sensitivity") +
      facet_wrap(. ~ method, ncol = 5) +
      theme(legend.position = "none") +
      theme(strip.text.x = element_text(size = 10))
    ggsave(file.path("output", "images", paste0("DE_p",p,"_corrp",corrp,"_all_models.png")),
           units = "in",
           height = 2.5,
           width = 10)

    if(p == 5000 & corrp == 0.5) {
      ggplot(plot_data, aes(x = fpr, y = tpr, color = log(delta_mean_v2))) +
        geom_point(size = 1, alpha = 0.5) +
        scale_color_gradient2(low = "blue", mid = "#BBBBBB", high = "red", midpoint = mean_log_diff) +
        xlim(c(0,0.75)) +
        ylim(c(0.25,1)) +
        xlab("1 - specificity") +
        ylab("sensitivity") +
        labs(color = "Log\nfold\nchange") +
        facet_wrap(. ~ method, ncol = 5) +
        theme(strip.text.x = element_text(size = 10))
      ggsave(file.path("output", "images", paste0("DE_p",p,"_corrp",corrp,"_all_models_LFC.png")),
             units = "in",
             height = 2.5,
             width = 10.5)
      
      
    }
  }
}

# ------------------------------------------------------------------------------
#   TPR x FPR with baseline + spike-in
# ------------------------------------------------------------------------------

p <- 100 # 5000
corrp <- 0 # 0.5
data <- readRDS(file.path("output", paste0("simresults_p",p,"_corrp",corrp,"_all.rds")))

plot_data <- data %>%
  select(delta_mean_v1, delta_mean_v2, rate, rate_type, method) %>%
  filter(method %in% c("baseline", "spike_in")) %>%
  pivot_wider(names_from = rate_type, values_from = rate)
# Assign factor order...
plot_data$method <- factor(plot_data$method, levels = c("baseline",
                                                        "spike_in"))
# ...then rename
levels(plot_data$method) <- c("NB GLM", "NB GLM + spike-in")

plot_data <- eliminate_outliers(plot_data)
log_diff <- plot_data$log_diff
mean_log_diff <- plot_data$mean_log_diff
plot_data <- plot_data$df

ggplot(plot_data, aes(x = fpr, y = tpr, color = log(delta_mean_v2))) +
  geom_point(size = 1, alpha = 0.5) +
  scale_color_gradient2(low = "blue", mid = "#BBBBBB", high = "red", midpoint = mean_log_diff) +
  xlim(c(0,0.75)) +
  ylim(c(0.25,1)) +
  xlab("1 - specificity") +
  ylab("sensitivity") +
  labs(color = "Log\nfold\nchange") +
  facet_wrap(. ~ method, ncol = 5) +
  theme(strip.text.x = element_text(size = 10))
ggsave(file.path("output", "images", paste0("DE_p",p,"_corrp",corrp,"_spikein_LFC.png")),
       units = "in",
       height = 2.5,
       width = 5.5)

# ------------------------------------------------------------------------------
#   TPR x FPR with baseline + spike-in (LOCALLY GENERATED DATA)
# ------------------------------------------------------------------------------

output_dir <- "p100_corrp0"
n <- 10
data <- readRDS(file.path("output", output_dir, paste0("simresults_p100_simulated_all.rds")))

plot_data <- data %>%
  select(uuid, delta_mean_v1, delta_mean_v2, rate, rate_type, method, partial_info, gap_totals) %>%
  pivot_wider(names_from = rate_type, values_from = rate)
# Assign factor order...
plot_data$method <- factor(plot_data$method, levels = c("baseline",
                                                        "DESeq2"))
# ...then rename
levels(plot_data$method) <- c("NB GLM", "DESeq2")

plot_data <- eliminate_outliers(plot_data)
log_diff <- plot_data$log_diff
mean_log_diff <- plot_data$mean_log_diff
plot_data <- plot_data$df

# Add deltas in totals for 2nd set of observed counts (the partial info ones!)
plot_data$delta_mean_v1_partial <- numeric(nrow(plot_data))
plot_data$delta_mean_v2_partial <- numeric(nrow(plot_data))
uuids <- plot_data %>% distinct(uuid) %>% pull(uuid)
for(uuid in uuids) {
  sim_data <- readRDS(file.path("output", output_dir, paste0(uuid, ".rds")))
  r1 <- rowSums(sim_data$observed_counts2[1:n,])
  r2 <- rowSums(sim_data$observed_counts2[(n+1):(n*2),])
  m1 <- mean(r1)
  m2 <- mean(r2)
  delta_mean_v1 <- max(c(m1, m2)) - min(c(m1, m2)) # absolute difference
  delta_mean_v2 <- max(c(m1, m2)) / min(c(m1, m2)) # absolute fold change
  plot_data[plot_data$uuid == uuid,]$delta_mean_v1_partial <- delta_mean_v1
  plot_data[plot_data$uuid == uuid,]$delta_mean_v2_partial <- delta_mean_v2
}

for(method in levels(plot_data$method)) {
  sub_plot_data <- plot_data[plot_data$method == method,]
  save_name <- paste0("performance_",
                      ifelse(method == "NB GLM", "NBGLM", method))
  
  partial_data <- sub_plot_data[sub_plot_data$partial_info == TRUE,]
  baseline_data <- sub_plot_data[sub_plot_data$partial_info == FALSE,]
  
  # Baseline data
  ggplot(baseline_data, aes(x = fpr, y = tpr, color = log(delta_mean_v2))) +
    geom_point(size = 2, alpha = 0.66) +
    scale_color_gradient2(low = "navy", mid = "#cccccc", high = "red",
                          midpoint = mean(log(baseline_data$delta_mean_v2))) +
    xlim(c(0,0.75)) +
    ylim(c(0.25,1)) +
    xlab("1 - specificity") +
    ylab("sensitivity") +
    labs(color = "Log\nfold\nchange") +
  ggsave(file.path("output", "images", paste0(save_name, "_baseline.png")),
         units = "in",
         height = 4,
         width = 5)
  # Partial info-augmented data
  partial_data$gap_restored <- partial_data$delta_mean_v2_partial / 
    partial_data$delta_mean_v2
  subset_partial_df <- partial_data[partial_data$gap_restored <= 1.0,]
  ggplot(subset_partial_df, aes(x = fpr, y = tpr, color = gap_restored)) +
    geom_point(size = 2, alpha = 0.66) +
    scale_color_gradient2(low = "#ff8c00", mid = "#cccccc", high = "#2abd4e",
                          midpoint = mean(subset_partial_df$gap_restored)) +
    xlim(c(0,0.75)) +
    ylim(c(0.25,1)) +
    xlab("1 - specificity") +
    ylab("sensitivity") +
    labs(color = "Percent\nfold change\nretained") +
    ggsave(file.path("output", "images", paste0(save_name, "_partial.png")),
           units = "in",
           height = 4,
           width = 5)
  
  # ----------------------------------------------------------------------------
  #   Plot baseline + partial into w/ discrete LOW/MED/HIGH LFC coloring
  # ----------------------------------------------------------------------------

  # n_levels <- 3 # cut into factor by tertiles
  # cuts <- quantile(sub_plot_data$delta_mean_v2,
  #                  probs = seq(from = 0, to = 1, length.out = n_levels+1))
  # cuts[1] <- -Inf
  # cuts[length(cuts)] <- Inf
  # delta_values <- sub_plot_data$delta_mean_v2
  # delta_values[sub_plot_data$partial_info == TRUE] <- sub_plot_data[sub_plot_data$partial_info == TRUE,]$delta_mean_v2_partial
  # lfc_factor <- cut(delta_values, breaks = cuts)
  # levels(lfc_factor) <- c("low", "med", "high")
  # palette_lfc <- colorRampPalette(c("navy", "#bbbbbb", "red"))(n_levels)
  # names(palette_lfc) <- levels(lfc_factor)
  # sub_plot_data$lfc_factor <- lfc_factor
  # 
  # partial_data <- sub_plot_data[sub_plot_data$partial_info == TRUE,]
  # baseline_data <- sub_plot_data[sub_plot_data$partial_info == FALSE,]
  # 
  # # Baseline data
  # ggplot(baseline_data, aes(x = fpr, y = tpr, color = lfc_factor)) +
  #   geom_point(size = 2, alpha = 0.66) +
  #   scale_color_manual(values = palette_lfc) +
  #   xlim(c(0,0.75)) +
  #   ylim(c(0.25,1)) +
  #   xlab("FPR") +
  #   ylab("TPR") +
  #   labs(title = paste0(method, " on observed counts")) +
  #   theme(legend.position = "none")
  # ggsave(file.path("output", "images", paste0(save_name, "_baseline.png")),
  #        units = "in",
  #        height = 4,
  #        width = 4)
  # # Partial info-augmented data
  # ggplot(partial_data, aes(x = fpr, y = tpr, color = lfc_factor)) +
  #   geom_point(size = 2, alpha = 0.66) +
  #   scale_color_manual(values = palette_lfc) +
  #   xlim(c(0,0.75)) +
  #   ylim(c(0.25,1)) +
  #   xlab("FPR") +
  #   ylab("TPR") +
  #   labs(color = "Log\nfold\nchange",
  #        title = paste0(method, " on partially informative counts")) +
  # ggsave(file.path("output", "images", paste0(save_name, "_partial.png")),
  #        units = "in",
  #        height = 4,
  #        width = 5)
  
  # ----------------------------------------------------------------------------
  #   Plot improvement in FPR with increasing gap closed
  # ----------------------------------------------------------------------------
  
  plot_subset <- sub_plot_data %>%
    filter(partial_info == TRUE)
  plot_df <- data.frame(percent_gap_restored = plot_subset$delta_mean_v2_partial /
                          plot_subset$delta_mean_v2,
                        fpr = plot_subset$fpr,
                        gap_orig = plot_subset$delta_mean_v2)
  ggplot(plot_df, aes(x = percent_gap_restored, y = fpr, color = log(gap_orig))) +
    geom_smooth(color = "#bbbbbb", alpha = 0.2) +
    geom_point() +
    scale_color_gradient2(low = "navy",
                          mid = "#bbbbbb",
                          high = "red",
                          midpoint = mean(log(plot_df$gap_orig))) +
    xlim(c(min(plot_df$percent_gap_restored), 1)) +
    labs(color = "Log\nfold\nchange",
         x = "Percent fold change retained",
         y = "1 - specificity")
  ggsave(file.path("output", "images", paste0(save_name, "_gap-vs-fpr.png")),
         units = "in",
         height = 4,
         width = 5)
}

# "Improvement" plots - TBD
# this_method <- "NB GLM"
this_method <- "DESeq2"
deltas <- c()
gaps <- c()
for(uuid in unique(plot_data$uuid)) {
  partial_entry <- plot_data[plot_data$uuid == uuid &
                               plot_data$method == this_method &
                               plot_data$partial_info == TRUE,]$fpr
  base_entry <- plot_data[plot_data$uuid == uuid &
                            plot_data$method == this_method &
                            plot_data$partial_info == FALSE,]$fpr
  deltas <- c(deltas, partial_entry - base_entry)
  gaps <- c(gaps, plot_data[plot_data$uuid == uuid &
                              plot_data$method == this_method &
                              plot_data$partial_info == TRUE,]$gap_totals)
}
plot_df <- data.frame(gaps = gaps, deltas = deltas)
ggplot(plot_df, aes(x = gaps, y = deltas)) +
  geom_point() +
  geom_smooth(method = "lm")

# ------------------------------------------------------------------------------
#   LFC distributions
# ------------------------------------------------------------------------------

lfc_data <- lfc_data %>%
  distinct(lfc, p, corrp)
lfc_data$corrp <- factor(lfc_data$corrp)
lfc_data$p <- factor(lfc_data$p)

# Plot fold change distribution: uncorrelated features, 0.1K/5K features
ggplot() +
  geom_histogram(data = lfc_data[lfc_data$p == 100 & lfc_data$corrp == 0,],
                 mapping = aes(x = lfc, fill = p),
                 bins = 20,
                 color = "white",
                 alpha = 0.5) +
  geom_histogram(data = lfc_data[lfc_data$p == 5000 & lfc_data$corrp == 0,],
                 mapping = aes(x = lfc, fill = p),
                 bins = 20,
                 color = "white",
                 alpha = 0.5) +
  scale_fill_manual(values = c("blue", "red")) +
  labs(x = "log fold change (total abundance)",
       fill = "Number of\nfeatures")
ggsave(file.path("output", "images", "LFC_histogram_100_vs_5000.png"),
       units = "in",
       height = 4,
       width = 5)

# Plot fold change distribution: 5K features, uncorrelated/correlated
ggplot() +
  geom_histogram(data = lfc_data[lfc_data$p == 5000 & lfc_data$corrp == 0,],
                 mapping = aes(x = lfc, fill = corrp),
                 bins = 20,
                 color = "white",
                 alpha = 0.5) +
  geom_histogram(data = lfc_data[lfc_data$p == 5000 & lfc_data$corrp == 0.5,],
                 mapping = aes(x = lfc, fill = corrp),
                 bins = 20,
                 color = "white",
                 alpha = 0.5) +
  scale_fill_manual(values = c("red", "blue")) +
  labs(x = "log fold change (total abundance)",
       fill = "% Correlated\nfeatures")
ggsave(file.path("output", "images", "LFC_histogram_uncorrelated_vs_correlated.png"),
       units = "in",
       height = 4,
       width = 5)

# ------------------------------------------------------------------------------
#   ROC curves with validation results superimposed
# ------------------------------------------------------------------------------

# Note: There are edgeR results in here but I think they are wrong. (They
# perfectly resemble the NB GLM results.)

dataset <- "Barlow" # Athanasiadou or Barlow

if(dataset == "Athanasiadou") {
  data <- readRDS(file.path("output", paste0("simresults_p5000_corrp0.5_all.rds")))
  method_list <- c("baseline", "DESeq2", "MAST", "ALDEx2")
  data <- data %>%
    filter(method %in% method_list) %>%
    filter(delta_mean_v2 < 2)
  method_labels <- c("NB GLM", "DESeq2", "MAST", "ALDEx2")
  validation <- readRDS(file.path("output", "Athanasiadou_validation_results.rds"))
  validation <- validation %>%
    filter(method %in% method_list)
  save_name <- "Athanasiadou_validation.png"
} else {
  data <- readRDS(file.path("output", paste0("simresults_p100_corrp0.5_all.rds")))
  method_list <- c("baseline", "DESeq2", "MAST", "ALDEx2")
  data <- data %>%
    filter(method %in% method_list) %>%
    filter(delta_mean_v2 > 3 & delta_mean_v2 < 4)
  method_labels <- c("NB GLM", "DESeq2", "MAST", "ALDEx2")
  validation <- readRDS(file.path("output", "Barlow_validation_results.rds"))
  validation <- validation %>%
    filter(method %in% method_list)
  save_name <- "Barlow_validation.png"
}

plot_data <- data %>%
  select(delta_mean_v2, rate, rate_type, method) %>%
  filter(method %in% method_list) %>%
  pivot_wider(names_from = rate_type, values_from = rate)
plot_data$method <- factor(plot_data$method)

plot_list <- list()
for(i in 1:length(method_list)) {
  method <- method_list[[i]]
  pl <- ggplot(plot_data[plot_data$method == method,], aes(x = fpr, y = tpr)) +
    geom_point(alpha = 0.5) +
    # geom_density_2d(color = "black", size = 0.5, alpha = 0.5) +
    geom_point(data = validation[validation$method == method,],
               aes(x = fpr, y = tpr), size = 3, color = "red") +
    xlim(c(0,0.75)) +
    ylim(c(0.25,1)) +
    labs(title = method_labels[i], x = "FPR", y = "TPR")
  plot_list[[i]] <- pl
}

pl <- grid.arrange(grobs = plot_list, ncol = 4)

ggsave(file.path("output", "images", save_name),
       plot = pl,
       units = "in",
       height = 2.75,
       width = 10)

# ------------------------------------------------------------------------------
#   Improvement due to partial total abundance info
# ------------------------------------------------------------------------------
# 
# p <- 100
# corrp_all <- c(0, 0.5)
# for(corrp in corrp_all) {
#   data <- readRDS(file.path("output", paste0("simresults_p",p,"_corrp",corrp,"_all.rds")))
#   data <- data %>%
#     filter(method %in% c("baseline", "spike_in")) %>%
#     filter(rate_type == "fpr") %>%
#     filter(delta_mean_v2 > 2)
#   
#   plot_data <- data.frame(cor_totals = c(), delta_fpr = c())
#   for(i in seq(from = 1, to = (nrow(data)-1), by = 2)) {
#     delta_fpr <- data$rate[i+1] - data$rate[i]
#     cor_totals <- data$cor_totals[i+1]
#     plot_data <- rbind(plot_data,
#                        data.frame(cor_totals = cor_totals, delta_fpr = delta_fpr))
#   }
#   plot_data <- plot_data %>%
#     filter(cor_totals >= 0)
#   ggplot(plot_data, aes(x = cor_totals, y = delta_fpr)) +
#     geom_point(size = 1) +
#     geom_smooth(method = "lm", color = "#888888", alpha = 0.3) +
#     labs(x = "correlation true/observed totals",
#          y = "change in FPR (w/ correlated totals)")
#   ggsave(file.path("output",
#                    "images",
#                    paste0("partial_total_performance_",p,"_corrp",corrp,".png")),
#          units = "in",
#          height = 3,
#          width = 3)
# }

# ------------------------------------------------------------------------------
#   ALDEx2 vs. DESeq2
# ------------------------------------------------------------------------------

# aldex2_data <- cbind(readRDS(file.path("output", paste0("simresults_p100_corrp0_all.rds"))),
#                      setting = "p100_corrp0")
# aldex2_data <- rbind(aldex2_data,
#                      cbind(readRDS(file.path("output", paste0("simresults_p5000_corrp0_all.rds"))),
#                            setting = "p5000_corrp0"))
# aldex2_data <- rbind(aldex2_data,
#                      cbind(readRDS(file.path("output", paste0("simresults_p5000_corrp0.5_all.rds"))),
#                            setting = "p5000_corrp0.5"))
# 
# aldex2_data <- aldex2_data %>%
#   filter(method == "ALDEx2") %>%
#   pivot_wider(names_from = rate_type, values_from = rate)
# 
# log_diff <- log(aldex2_data$delta_mean_v2)
# mean_log_diff <- mean(log_diff)
# limit <- 3*sd(log_diff)
# idx_retain <- which(log_diff > mean_log_diff - limit & log_diff < mean_log_diff + limit)
# aldex2_data <- aldex2_data[idx_retain,]
# 
# # unary_palette <- palette[6]
# 
# case.labs <- c("ALDEx2, 100 uncorrelated feat.", "ALDEx2, 5000 uncorrelated feat.", "ALDEx2, 5000 correlated feat.")
# names(case.labs) <- c("p100_corrp0", "p5000_corrp0", "p5000_corrp0.5")
# 
# ggplot(aldex2_data, aes(x = fpr, y = tpr, color = log(delta_mean_v2))) +
#   geom_point(size = 2, alpha = 0.5) +
#   scale_color_gradient2(low = "blue", mid = "#BBBBBB", high = "red", midpoint = mean_log_diff) +
#   xlim(c(0,0.75)) +
#   ylim(c(0,1)) +
#   xlab("FPR") +
#   ylab("TPR") +
#   facet_wrap(. ~ setting, labeller = labeller(setting = case.labs)) +
#   theme(legend.position = "none")
# ggsave(file.path("output",
#                  "images",
#                  "ALDEx2_performance_LFC.png"),
#        units = "in",
#        height = 3,
#        width = 7)

# ------------------------------------------------------------------------------
#   DESeq2 performance (FPR) as a function of log fold difference in means
# ------------------------------------------------------------------------------

# lfc_plot <- TRUE
# 
# if(lfc_plot) {
#   unary_palette <- palette[2]
#   
#   pp <- c(1000, 5000, 15000)
#   data <- cbind(p = 1000, readRDS(file.path("output", paste0("simresults_p",pp[1],"_simulated_all.rds"))))
#   for(p in pp[2:length(pp)]) {
#     data <- rbind(data, cbind(p = p, readRDS(file.path("output", paste0("simresults_p",p,"_simulated_all.rds")))))
#   }
#   data <- data %>%
#     filter(method == "DESeq2") %>%
#     filter(p > 100) %>%
#     filter(rate_type == "fpr") %>%
#     filter(delta_mean_v2 <= 10) %>%
#     select(p, delta_mean_v2, rate, rate_type, method)
#   data$p <- factor(data$p)
#   
#   ggplot(data, aes(x = log(delta_mean_v2), y = rate, color = rate_type)) + # plot an ROC curve
#     geom_point(size = 1, alpha = 0.25) +
#     geom_smooth(method = "loess", alpha = 0.5) +
#     scale_color_manual(values = unary_palette) +
#     theme(legend.position = "none") +
#     ylim(c(0,0.15))
#   ggsave(file.path("output", "images", "DESeq2_performance_curves.png"),
#          units = "in",
#          height = 3,
#          width = 3)
# }

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



