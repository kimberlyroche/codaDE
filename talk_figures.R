library(tidyverse)
library(cowplot)
library(ggridges)
library(ggrepel)

output_dir <- "C:/Users/kimbe/Desktop"

mpalette <- c("#46A06B", "#FF5733", "#8065f0")
# names(mpalette) <- c("scran", "DESeq2", "ALDEx2")
names(mpalette) <- c("ALDEx2", "DESeq2", "scran")

dpalette <- colorRampPalette(brewer.pal(12, "Paired"))(12)

# ------------------------------------------------------------------------------
#   Simulated toy count data
# ------------------------------------------------------------------------------

# Counts
dispersion <- 1000
plot_df <- data.frame(x = 1:10,
                      y = c(rnbinom(5, mu = 100, size = dispersion), rnbinom(5, mu = 50, size = dispersion)),
                      type = rep(c("control", "treatment"), each = 5))

# hcolor <- "#9b445d"
hcolor <- "#e38552"

p <- ggplot(plot_df, aes(x = factor(x), y = y, fill = type, color = type)) +
  geom_point(size = 5, shape = 21) +
  theme_classic() +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # ylim(c(50, 250)) +
  labs(fill = "",
       x = "sample",
       y = "abundance") +
  scale_fill_manual(values = rep(hcolor, 2)) +
  scale_color_manual(values = rep(hcolor, 2))

ggsave(file.path(output_dir, "counts.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 4,
       width = 4)

# Absolute and relative abundances

ref_obj <- build_simulated_reference(p = 100, log_mean = 1, log_var = 4,
                                     perturb_size = 1.5, base_correlation = NULL,
                                     concentration = 10000, save_name = NULL)
sim_obj <- simulate_sequence_counts(n = 5,
                                     p = 20,
                                     # ref_data_name = "simulated",
                                     data_obj = ref_obj,
                                     replicate_noise = 0,
                                     proportion_da = 0.75)

palette <- generate_highcontrast_palette(20)

plot_stacked_bars(sim_obj$abundances[c(2:5,1,6:10),], palette = palette)
ggsave(file.path(output_dir, "abs_ab.svg"),
       dpi = 100,
       units = "in",
       height = 4,
       width = 5)

cpm <- t(sim_obj$observed_counts2[c(2:5,1,6:10),])
cpm <- apply(cpm, 2, function(x) x/sum(x))
cpm <- t(cpm)
plot_stacked_bars(cpm, palette = palette)
ggsave(file.path(output_dir, "rel_ab.svg"),
       dpi = 100,
       units = "in",
       height = 4,
       width = 5)

# ------------------------------------------------------------------------------
#   ROC
# ------------------------------------------------------------------------------

p <- ggplot(data, aes(x = 1-FPR, y = TPR, fill = factor(P))) +
  geom_point(size = 2, shape = 21) +
  # geom_smooth(method = "loess") +
  theme_bw() +
  labs(x = "false positive rate",
       y = "true positive rate",
       fill = "feature number") +
  scale_fill_brewer(palette = "Reds") +
  theme(legend.position = "bottom")
pm <- ggMarginal(p, margins = "both", type = "histogram", bins = 50,
                 fill = "#888888", color = "#333333", size = 20)
ggsave(file.path(output_dir, "roc.svg"),
       pm,
       dpi = 100,
       units = "in",
       height = 6,
       width = 6)

# ------------------------------------------------------------------------------
#   ROC labeled for percent DE
# ------------------------------------------------------------------------------

temp <- data %>%
  left_join(data3, by = "UUID")

p <- ggplot(temp, aes(x = 1-FPR.x, y = TPR.x, fill = PERCENT_DIFF_REALIZ)) +
  geom_point(size = 2, shape = 21) +
  # geom_smooth(method = "loess") +
  theme_bw() +
  labs(x = "false positive rate",
       y = "true positive rate",
       fill = "percent differentially\nabundant features") +
  # scale_fill_distiller(palette = "Spectral") +
  scale_fill_gradient2(low = "#18BD6F", mid = "white", high = "#FF5733", midpoint = 0.5) +
  theme(legend.position = "bottom")
ggsave(file.path(output_dir, "roc.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 6,
       width = 6)

# ------------------------------------------------------------------------------
#   R2 for most predictive feature (FPR)
# ------------------------------------------------------------------------------

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
res <- dbGetQuery(conn, paste0("SELECT UUID, FW_CLR_PFC05_D FROM characteristics WHERE TYPE='cpm' AND partial=0;"))
dbDisconnect(conn)

data2 <- data %>%
  left_join(res, by = "UUID")

p <- ggplot(data2, aes(x = FW_CLR_PFC05_D, y = PERCENT_DIFF_REALIZ, fill = factor(P))) +
  geom_point(size = 2, shape = 21) +
  theme_bw() +
  labs(x = "percent of features with < 50% change in mean-relative log abundance",
       y = "percent of features differentially abundant",
       fill = "feature number") +
  scale_fill_brewer(palette = "Reds") +
  theme(legend.position = "bottom")
ggsave(file.path(output_dir, "pred_feat.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 6,
       width = 6)

# ------------------------------------------------------------------------------
#   Identifiable "problems" in simulated data
# ------------------------------------------------------------------------------

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
char_data <- dbGetQuery(conn, paste0("SELECT * FROM characteristics WHERE TYPE='cpm' AND partial=0;"))
dbDisconnect(conn)

data3 <- char_data %>%
  left_join(data2 %>% select(UUID, P, METHOD, FPR, TPR), by = "UUID") %>%
  filter(!is.na(P))

# This takes ~17 min.
start <- Sys.time()
data3$predicted <- NA
for(i in 1:nrow(data3)) {
  if(i %% 100 == 0) {
    cat(paste0(i, " / ", nrow(data3), "\n"))
  }
  if(data3$METHOD[i] == "ALDEx2") {
    data3$predicted[i] <- predict(fit_a$result, data3[i,])
  }
  if(data3$METHOD[i] == "DESeq2") {
    data3$predicted[i] <- predict(fit_d$result, data3[i,])
  }
  if(data3$METHOD[i] == "scran") {
    data3$predicted[i] <- predict(fit_s$result, data3[i,])
  }
}
end <- Sys.time()
diff <- end - start
print(diff)

problems <- rep(FALSE, 1:nrow(use_data))

p <- ggplot() +
  geom_point(data = data3 %>% filter(predicted >= 0.9),
             mapping = aes(x = 1-FPR, y = TPR),
             color = "#dddddd") +
  geom_point(data = data3 %>% filter(predicted < 0.9),
             mapping = aes(x = 1-FPR, y = TPR),
             size = 2, shape = 21, fill = "#4CB1EC", alpha = 0.5) +
  geom_segment(data = data.frame(x = 0.1, xend = 0.1, y = 0, yend = 1),
               mapping = aes(x = x, xend = xend, y = y, yend = yend),
               size = 1,
               linetype = "dashed") +
  # geom_smooth(method = "loess") +
  theme_bw() +
  labs(x = "false positive rate",
       y = "true positive rate") +
  theme(legend.position = "bottom")
ggsave(file.path(output_dir, "pred.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 6,
       width = 6)

# What proportion of data sets are high FPR
n <- data3 %>% filter(predicted >= 0.9) %>% tally() %>% pull(n)
n/nrow(data3)

# ------------------------------------------------------------------------------
#   R2 for prediction
# ------------------------------------------------------------------------------

fit_a <- readRDS("output/predictive_fits/oracle/regression_cpm/FPR_ALDEx2.rds")
plot_df <- data.frame(method = "ALDEx2",
                      response = fit_a$test_response,
                      prediction = predict(fit_a$result, fit_a$test_features))

fit_d <- readRDS("output/predictive_fits/oracle/regression_cpm/FPR_DESeq2.rds")
plot_df <- rbind(plot_df,
                 data.frame(method = "DESeq2",
                            response = fit_d$test_response,
                            prediction = predict(fit_d$result, fit_d$test_features)))

fit_s <- readRDS("output/predictive_fits/oracle/regression_cpm/FPR_scran.rds")
plot_df <- rbind(plot_df,
                 data.frame(method = "scran",
                            response = fit_s$test_response,
                            prediction = predict(fit_s$result, fit_s$test_features)))

p <- ggplot(plot_df, aes(x = 1-response, y = 1-prediction, fill = data3)) +
  geom_point(size = 2, shape = 21) +
  theme_bw() +
  labs(x = "false positive rate (observed)",
       y = "false positive rate (predicted)") +
  scale_fill_manual(values = mpalette) +
  theme(legend.position = "bottom") +
  xlim(c(0,1)) +
  ylim(c(0,1))
ggsave(file.path(output_dir, "pred.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 6.5,
       width = 6)

# ------------------------------------------------------------------------------
#   Histograms
# ------------------------------------------------------------------------------

# First load `data` object from `figure_1-2-S2-S3.R`
data$METHOD <- factor(data$METHOD, levels = c("scran", "DESeq2", "ALDEx2"))

# Set up a factor for percent of moving features
data <- data %>%
  mutate(PDIFF_DISCRETE = case_when(
    PERCENT_DIFF_REALIZ < 0.1 ~ "0-10%",
    PERCENT_DIFF_REALIZ < 0.25 ~ "10-25%",
    PERCENT_DIFF_REALIZ < 0.5 ~ "25-50%",
    TRUE ~ ">50%"
  ))
data$PDIFF_DISCRETE <- factor(data$PDIFF_DISCRETE,
                              levels = c("0-10%","10-25%", "25-50%", ">50%"))

# Methods are similar
# Long tail of "bad" outcomes
ggplot(data, aes(x = 1-FPR, y = METHOD, fill = METHOD)) +
  geom_density_ridges(stat = "binline", bins = 40, scale = 1) +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  scale_y_discrete(expand = expansion(add = c(0.2, 1.1))) +
  coord_cartesian(clip = "off") +
  guides(fill = "none") +
  labs(x = "false positive rate") +
  scale_fill_manual(values = mpalette)
ggsave(file.path(output_dir, "temp.svg"),
       dpi = 100,
       units = "in",
       height = 4,
       width = 4)

# Role of P
ggplot(data, aes(x = factor(P), y = 1-FPR)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_bw() +
  labs(x = "feature number",
       y = "false positive rate")

plots <- list()
no_y <- FALSE
no_x <- TRUE
for(this_method in c("ALDEx2", "DESeq2", "scran")) {
  # Role of percent differential features
  p <- ggplot(data %>% filter(METHOD == this_method), aes(x = 1-FPR, y = PDIFF_DISCRETE)) +
    geom_density_ridges(stat = "binline", bins = 40, scale = 0.95) +
    theme_bw() +
    scale_y_discrete(expand = expansion(add = c(0.15, 0.4))) +
    coord_cartesian(clip = "off") +
    labs(x = ifelse(no_x, "\n", "\nfalse positive rate"),
         y = ifelse(no_y, "", "percent differentially abundant features\n"),
         title = this_method)
  if(no_y) {
    p <- p +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
  plots[[length(plots)+1]] <- p
  no_x <- !no_x
  no_y <- TRUE
}
plot_grid(plotlist = plots, ncol = 3)

# Role of FC
plots <- list()
no_y <- FALSE
no_x <- TRUE
for(this_method in c("ALDEx2", "DESeq2", "scran")) {
  # Role of percent differential features
  p <- ggplot(data %>% filter(METHOD == this_method), aes(x = 1-FPR, y = FC_plot)) +
    geom_density_ridges(stat = "binline", bins = 40, scale = 0.95) +
    theme_bw() +
    scale_y_discrete(expand = expansion(add = c(0.15, 0.4))) +
    coord_cartesian(clip = "off") +
    labs(x = ifelse(no_x, "\n", "\nfalse positive rate"),
         y = ifelse(no_y, "", "percent differentially abundant features\n"),
         title = this_method)
  if(no_y) {
    p <- p +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
  plots[[length(plots)+1]] <- p
  no_x <- !no_x
  no_y <- TRUE
}
plot_grid(plotlist = plots, ncol = 3, rel_widths = c(1.5, 1, 1))

# ------------------------------------------------------------------------------
#   Real data calls superimposed on simulated background
# ------------------------------------------------------------------------------

real_calls <- NULL
real_files <- list.files(path = "output/real_data_calls/no_norm", pattern = ".*\\.rds", full.names = TRUE)
for(file in real_files) {
  calls <- readRDS(file)
  pieces <- str_split(str_split(str_split(file, "\\/")[[1]][4], "\\.")[[1]][1], "_")[[1]]
  real_calls <- rbind(real_calls,
                      data.frame(method = pieces[3],
                                 dataset = pieces[4],
                                 tpr = calls$rates$TPR,
                                 fpr = calls$rates$FPR,
                                 fp = sum(calls$rates$FP_calls),
                                 tp = sum(calls$rates$TP_calls),
                                 true_diff = sum(p.adjust(calls$all_calls$oracle_calls, method = "BH") < 0.05),
                                 n = length(calls$all_calls$oracle_calls)))
}

temp <- real_calls %>%
  mutate(true_diff_percent = true_diff/n) %>%
  arrange(true_diff_percent) %>%
  mutate(index = 1:n())

temp <- temp %>%
  mutate(datatype = case_when(
    dataset %in% c("Barlow", "VieiraSilva") ~ "16S",
    dataset %in% c("Hagai", "Monaco", "Song") ~ "bulk",
    TRUE ~ "sc"
  ))

p0 <- ggplot(data.frame(x = c(1,2,3), y = c(1,2,3), method = c("ALDEx2", "DESeq2", "scran")),
             aes(x = x, y = y, shape = method)) +
  geom_point(size = 3) +
  theme(legend.position = "bottom")
shape_legend <- get_legend(p0)

point_sz <- 4
jitter <- 0
p <- ggplot() +
  geom_jitter(data = temp %>% filter(method == "ALDEx2"),
              mapping = aes(x = true_diff_percent*100, y = fp/n*100, fill = dataset),
              size = point_sz,
              shape = 21,
              width = jitter) +
  scale_fill_manual(values = dpalette) +
  theme(legend.position = "bottom")
legend <- get_legend(p)
p <- p +
  geom_jitter(data = temp %>% filter(method == "DESeq2"),
              mapping = aes(x = true_diff_percent*100, y = fp/n*100, fill = dataset),
              size = point_sz - 0.5,
              shape = 24,
              width = jitter) +
  theme(legend.position = "none")
p <- p +
  geom_jitter(data = temp %>% filter(method == "scran"),
              mapping = aes(x = true_diff_percent*100, y = fp/n*100, fill = dataset),
              size = point_sz,
              shape = 22,
              width = jitter) +
  theme_bw() +
  labs(x = "percent differential features",
       y = "false positives as percent total features") +
  theme(legend.position = "none") +
  xlim(c(0, 100))

p_allreal <- plot_grid(p, legend, shape_legend, ncol = 1, rel_heights = c(1, 0.2, 0.1))

ggsave(file.path(output_dir, "allreal.svg"),
       p_allreal,
       dpi = 100,
       units = "in",
       height = 7,
       width = 7)

type <- "sc"
temp$alpha <- 1
temp$alpha[temp$datatype != type] <- 0.5
p <- ggplot() +
  geom_jitter(data = temp %>% filter(method == "ALDEx2"),
              mapping = aes(x = true_diff_percent*100, y = fp/n*100, fill = dataset, alpha = alpha),
              size = point_sz,
              shape = 21,
              width = jitter) +
  scale_fill_manual(values = dpalette) +
  scale_alpha(guide = 'none') +
  theme(legend.position = "bottom")
legend <- get_legend(p)
p <- p +
  geom_jitter(data = temp %>% filter(method == "DESeq2"),
              mapping = aes(x = true_diff_percent*100, y = fp/n*100, fill = dataset, alpha = alpha),
              size = point_sz - 0.5,
              shape = 24,
              width = jitter) +
  theme(legend.position = "none")
p <- p +
  geom_jitter(data = temp %>% filter(method == "scran"),
              mapping = aes(x = true_diff_percent*100, y = fp/n*100, fill = dataset, alpha = alpha),
              size = point_sz,
              shape = 22,
              width = jitter) +
  theme_bw() +
  labs(x = "percent differential features",
       y = "false positives as percent total features") +
  theme(legend.position = "none") +
  xlim(c(0, 100))

p_allreal <- plot_grid(p, legend, shape_legend, ncol = 1, rel_heights = c(1, 0.2, 0.1))

ggsave(file.path(output_dir, paste0("allreal_", type, ".svg")),
       p_allreal,
       dpi = 100,
       units = "in",
       height = 7,
       width = 7)





p <- ggplot() +
  geom_point(data = data3,
             mapping = aes(x = 1-FPR, y = TPR),
             color = "#dddddd") +
  geom_point(data = real_calls,
             mapping = aes(x = fpr, y = tpr, fill = method),
             size = 3, shape = 21, fill = "#4CB1EC") +
  # geom_text_repel(data = real_calls,
  #                 mapping = aes(x = fpr, y = tpr, label = dataset),
  #                 size = 4) +
  theme_bw() +
  labs(x = "false positive rate",
       y = "true positive rate") +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "real.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 6,
       width = 6)

p0 <- ggplot(data.frame(x = c(1,2,3), y = c(1,2,3), method = c("ALDEx2", "DESeq2", "scran")),
             aes(x = x, y = y, shape = method)) +
  geom_point(size = 3) +
  theme(legend.position = "bottom")
shape_legend <- get_legend(p0)

splits <- list(
  list(datasets = c("Hashimshony", "Monaco", "Kimmerling"),
       indices = c(1:3)),
  list(datasets = c("VieiraSilva", "Muraro", "Hagai", "Gruen"),
       indices = c(c(4,5,10,7))),
  list(datasets = c("Song", "Barlow", "Yu"),
       indices = c(c(7,8,9))),
  list(datasets = c("Klein", "Owens"),
       indices = c(11:12))
)

for(i in 1:length(splits)) {
  p <- ggplot() +
    geom_point(data = data3,
               mapping = aes(x = 1-FPR, y = TPR),
               color = "#dddddd") +
    geom_point(data = real_calls %>% filter(dataset %in% splits[[i]]$datasets & method == "ALDEx2"),
               mapping = aes(x = fpr, y = tpr, fill = dataset),
               size = 3, shape = 21, stroke = 1) +
    scale_fill_manual(values = dpalette[splits[[i]]$indices]) +
    theme(legend.position = "bottom")
  fill_legend <- get_legend(p)
  
  p <- p +
    geom_point(data = real_calls %>% filter(dataset %in% splits[[i]]$datasets & method == "DESeq2"),
               mapping = aes(x = fpr, y = tpr, fill = dataset),
               size = 3, shape = 24, stroke = 1) +
    geom_point(data = real_calls %>% filter(dataset %in% splits[[i]]$datasets & method == "scran"),
               mapping = aes(x = fpr, y = tpr, fill = dataset),
               size = 3, shape = 22, stroke = 1) +
    theme_bw() +
    labs(x = "false positive rate",
         y = "true positive rate") +
    theme(legend.position = "none")
  
  p_combo <- plot_grid(p, fill_legend, shape_legend, ncol = 1,
                       rel_heights = c(1, 0.08, 0.06))
  
  ggsave(file.path(output_dir, paste0("real", i, ".svg")),
         p_combo,
         dpi = 100,
         units = "in",
         height = 6.4,
         width = 6)
}

# ------------------------------------------------------------------------------
#   Real data calls superimposed on predicted intervals for FPR only
# ------------------------------------------------------------------------------

# First load `results` object from `figure_5-6.R`

results$DE_method <- factor(results$DE_method)
results$lower90[!is.na(results$lower90)] <- 1 - results$lower90[!is.na(results$lower90)]
results$lower50[!is.na(results$lower50)] <- 1 - results$lower50[!is.na(results$lower50)]
results$upper50[!is.na(results$upper50)] <- 1 - results$upper50[!is.na(results$upper50)]
results$upper90[!is.na(results$upper90)] <- 1 - results$upper90[!is.na(results$upper90)]
results$point[!is.na(results$point)] <- 1 - results$point[!is.na(results$point)]

# Replace the mean with the median for this data point.
# I should be using the median for all these boxplots but I goofed and recorded
#   the means. For most dataset x method combos, these are very similar but for
#   VS x scran, a few outliers drag the mean outside the IQR, making for a crazy
#   boxplot.
results[results$dataset == "VieiraSilva" &
          results$result_type == "predicted" &
          results$DE_method == "scran" &
          results$score_type == "FPR",]$point <- 1 - 0.9747143

plot_intervals <- function(this_dataset, results, label = NULL) {
  legend <- NULL
  # Wrangle the poorly organized data
  plot_df <- results %>%
    filter(dataset == this_dataset) %>%
    filter(threshold == 1) %>%
    filter(result_type == "predicted") %>%
    filter(score_type == "FPR") %>%
    dplyr::select(!c(threshold, score_type))
  plot_df$true <- NA
  for(j in 1:nrow(plot_df)) {
    plot_df$true[j] <- results %>%
      filter(dataset == this_dataset) %>%
      filter(threshold == 1) %>%
      filter(result_type == "true") %>%
      filter(score_type == "FPR") %>%
      filter(DE_method == plot_df$DE_method[j]) %>%
      pull(point)
    
    p <- ggplot() +
      geom_boxplot(data = plot_df,
                   mapping = aes(x = DE_method,
                                 # ymin = lower90,
                                 ymax = lower90,
                                 # lower = lower50,
                                 upper = lower50,
                                 middle = point,
                                 # upper = upper50,
                                 lower = upper50,
                                 # ymax = upper90,
                                 ymin = upper90,
                                 fill = DE_method),
                   stat = "identity", color = "#666666", width = 0.5, alpha = 0.4) +
      geom_point(data = plot_df,
                 mapping = aes(x = DE_method, y = true, fill = DE_method),
                 size = 4, shape = 21, stroke = 1) +
      ylim(c(0,1)) +
      scale_fill_manual(values = mpalette) +
      # scale_fill_brewer(palette = "Set2") +
      theme_bw() +
      theme(legend.position = "none",
            plot.title = element_text(size = 10)) +
      labs(x = "",
           y = "false positive rate",
           title = ifelse(is.null(label), paste0(this_dataset, " et al."), label))
  }
  p
}

plotlist <- list(Gruen = plot_intervals("Gruen", results, "Gruen et al. (single-cell)"),
                 Hagai = plot_intervals("Hagai", results, "Hagai et al. (bulk RNA-seq)"),
                 Muraro = plot_intervals("Muraro", results, "Muraro et al. (single-cell)"),
                 VieiraSilva = plot_intervals("VieiraSilva", results, "Vieira-Silva et al. (16S metabarcoding)"),
                 Klein = plot_intervals("Klein", results, "Klein et al. (single-cell)"),
                 Owens = plot_intervals("Owens", results, "Owens et al. (single-cell)"))

prow1 <- plot_grid(plotlist[[1]], plotlist[[2]], ncol = 2,
                   rel_widths = c(1, 1))
prow2 <- plot_grid(plotlist[[3]], plotlist[[4]], ncol = 2,
                   rel_widths = c(1, 1))
pleft <- plot_grid(prow1, prow2, ncol = 1)

pright <- plot_grid(plotlist[[5]], plotlist[[6]], ncol = 1)

ggsave(file.path(output_dir, "intervals.svg"),
       pleft,
       dpi = 100,
       units = "in",
       height = 6,
       width = 6)

ggsave(file.path(output_dir, "intervals2.svg"),
       pright,
       dpi = 100,
       units = "in",
       height = 6,
       width = 3)

# ------------------------------------------------------------------------------
#   False discovery rates for <25% DE 
# ------------------------------------------------------------------------------

data <- data %>%
  mutate(PDIFF_DISCRETE = case_when(
    PERCENT_DIFF_REALIZ < 0.1 ~ "0-10%",
    PERCENT_DIFF_REALIZ < 0.25 ~ "10-25%",
    PERCENT_DIFF_REALIZ < 0.5 ~ "25-50%",
    TRUE ~ ">50%"
  ))
data$PDIFF_DISCRETE <- factor(data$PDIFF_DISCRETE,
                              levels = c("0-10%","10-25%", "25-50%", ">50%"))
# fdrs <- numeric(nrow(data))
fprs <- numeric(nrow(data))
fake_yes <- numeric(nrow(data))
for(i in 1:nrow(data)) {
  if(i %% 1000 == 0) {
    cat(paste0("Iteration ", i, "\n"))
  }
  x <- as.numeric(str_split(data[i,]$ORACLE_BASELINE, ";")[[1]])
  y <- as.numeric(str_split(data[i,]$CALLS, ";")[[1]])
  
  x <- p.adjust(x, method = "BH")
  y <- p.adjust(y, method = "BH")
  
  x_d <- factor(x < 0.001)
  levels(x_d) <- c("not DE", "DE")
  y_d <- factor(y < 0.001)
  levels(y_d) <- c("not DE", "DE")
  
  # fdrs[i] <- sum(y_d == "DE" & x_d == "not DE")/(sum(y_d == "DE"))
  fprs[i] <- sum(y_d == "DE" & x_d == "not DE")/(sum(x_d == "not DE"))
  fake_yes[i] <- sum(y_d == "DE" & x_d == "not DE")
}
data$FPR1 <- fprs

data$FDR1[is.nan(data$FDR1)] <- 0

p <- ggplot(data, aes(x = FDR1, y = PDIFF_DISCRETE)) +
  geom_density_ridges(stat = "binline", bins = 40, scale = 0.95) +
  facet_wrap(. ~ METHOD) +
  theme_bw() +
  scale_y_discrete(expand = expansion(add = c(0.15, 0.4))) +
  coord_cartesian(clip = "off") +
  labs(x = "\nfalse positive rate",
       y = "percent differentially abundant features\n")

ggsave(file.path(output_dir, "fdr1.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 4,
       width = 8)

data$FP_n <- fake_yes

temp <- data
temp$P <- paste0(temp$P, " features")
temp$METHOD <- factor(temp$METHOD, levels = c("ALDEx2", "DESeq2", "scran"))
temp <- temp %>%
  mutate(PDIFF_DISCRETE = case_when(
    PERCENT_DIFF_REALIZ < 0.2 ~ "0-20%",
    PERCENT_DIFF_REALIZ < 0.4 ~ "20-40%",
    PERCENT_DIFF_REALIZ < 0.6 ~ "40-60%",
    PERCENT_DIFF_REALIZ < 0.8 ~ "60-80%",
    TRUE ~ "80-100%"
  ))
temp$PDIFF_DISCRETE <- factor(temp$PDIFF_DISCRETE,
                              levels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%"))

# p <- ggplot(temp, aes(x = PERCENT_DIFF_REALIZ, y = FP_n, fill = METHOD)) +
#   geom_point(size = 1.5, shape = 21) +
#   facet_wrap(. ~ METHOD + P, scales = "free_y") +
#   theme_bw() +
#   scale_fill_manual(values = mpalette) +
#   theme(legend.position = "none") +
#   labs(x = "\npercent differentially abundant features",
#        y = "false positive count\n")

temp$title <- paste0(temp$METHOD, " x ", temp$P)
p1 <- ggplot(temp %>% filter(P == "100 features"), aes(x = PDIFF_DISCRETE, y = FP_n, fill = METHOD)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(. ~ title, ncol = 1) + #, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = mpalette) +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  labs(x = "\npercent differentially abundant features",
       y = "false positive count\n") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p2 <- ggplot(temp %>% filter(P == "1000 features"), aes(x = PDIFF_DISCRETE, y = FP_n, fill = METHOD)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(. ~ title, ncol = 1) + #, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = mpalette) +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  labs(x = "\npercent differentially abundant features",
       y = "false positive count\n") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p3 <- ggplot(temp %>% filter(P == "5000 features"), aes(x = PDIFF_DISCRETE, y = FP_n, fill = METHOD)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(. ~ title, ncol = 1) + #, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = mpalette) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(x = "\npercent differentially abundant features",
       y = "false positive count\n") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p1_padded <- plot_grid(p1, NULL, ncol = 1, rel_heights = c(1, 0.05))
p3_padded <- plot_grid(p3, NULL, ncol = 1, rel_heights = c(1, 0.05))
p <- plot_grid(p1_padded, p2, p3_padded, ncol = 3, rel_widths = c(1, 0.85, 0.85))

ggsave(file.path(output_dir, "abs.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 6.5,
       width = 8)











