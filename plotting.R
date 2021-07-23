source("path_fix.R")

library(tidyverse)
library(codaDE)
library(RSQLite)
library(gridExtra)
library(RColorBrewer)
library(randomForest)

source("ggplot_fix.R")

label_combo <- function(data, plot_tag, logical_vec) {
  data$flag <- FALSE
  data$flag[logical_vec] <- TRUE
  plot_ROC(data, plot_tag, "flag", "selected")
}

plot_ROC <- function(data, plot_tag = NULL, fill_var = NULL, fill_var_label = NULL) {
  if(is.null(fill_var)) {
    fill_var <- "METHOD"
    fill_var_label <- "Method"
  }
  data$FPR <- 1 - data$FPR
  if(fill_var == "flag") {
    # Layer the points
    pl <- ggplot() +
      geom_point(data = data[data$flag == FALSE,],
                 mapping = aes_string(x = "FPR", y = "TPR"),
                 size = 2,
                 shape = 21,
                 fill = "#dddddd") +
      geom_point(data = data[data$flag == TRUE,],
                 mapping = aes_string(x = "FPR", y = "TPR"),
                 size = 2,
                 shape = 21,
                 fill = "#1ab079") +
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      labs(x = "specificity (1 - FPR)",
           y = "sensitivity (TPR)",
           fill = fill_var_label) +
      theme_bw() +
      facet_wrap(. ~ METHOD, ncol = 4)
  } else {
    pl <- ggplot(data, aes_string(x = "FPR", y = "TPR", fill = fill_var)) +
      geom_point(size = 2, shape = 21) +
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      labs(x = "specificity (1 - FPR)",
           y = "sensitivity (TPR)",
           fill = fill_var_label) +
      theme_bw() +
      facet_wrap(. ~ METHOD, ncol = 4)
  }
  if(is.numeric(data[[fill_var]])) {
    pl <- pl +
      scale_fill_distiller(palette = "RdYlBu")
  }
  if(is.factor(data[[fill_var]])) {
    pl <- pl +
      scale_fill_brewer(palette = "RdYlBu")
  }
  if(fill_var == "METHOD") {
    pl <- pl +
      scale_fill_manual(values = palette)
  }
  if(!is.null(plot_tag)) {
    ggsave(file.path("output", "images", paste0(fill_var,
                                                "_",
                                                plot_tag,
                                                ".png")),
           plot = pl,
           dpi = 100,
           units = "in",
           height = 3,
           width = 10)
  }
  show(pl)
}

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

palette <- c("#46A06B", "#FF5733", "#EF82BB", "#7E54DE", "#E3C012", "#B95D6E")

# ------------------------------------------------------------------------------
#   TPR x FPR plots
# ------------------------------------------------------------------------------

# Create directories manually
dir.create("output", showWarnings = FALSE)
dir.create(file.path("output", "images"), showWarnings = FALSE)
# dir.create(file.path("output", "images", "full_results"), showWarnings = FALSE)
# dir.create(file.path("output", "images", "full_results", "self_ref"), showWarnings = FALSE)
# dir.create(file.path("output", "images", "full_results", "oracle_ref"), showWarnings = FALSE)
# dir.create(file.path("output", "images", "full_results", "self_ref", "partial"), showWarnings = FALSE)
# dir.create(file.path("output", "images", "full_results", "self_ref", "no_partial"), showWarnings = FALSE)
# dir.create(file.path("output", "images", "full_results", "oracle_ref", "partial"), showWarnings = FALSE)
# dir.create(file.path("output", "images", "full_results", "oracle_ref", "no_partial"), showWarnings = FALSE)

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

ps <- c(100, 1000, 5000)
partials <- c(0, 1)
references <- c("oracle", "self")

qq_res <- dbGetQuery(conn, "SELECT FC_ABSOLUTE FROM datasets")$FC_ABSOLUTE
qq_res <- sapply(qq_res, function(x) {
  if(x < 1) {
    1 / x
  } else {
    x
  }
})
qq <- quantile(qq_res, probs = seq(from = 0, to = 1, length.out = 4))

# Build a fold change color scale
# ramped_palette1 <- colorRampPalette(c("#00cc00", "#ffffff", "#fe8b00"))(7)
# ramped_palette2 <- colorRampPalette(c("#7479c4", "#e63030"))(3)

for(p in ps) {
  for(partial in partials) {
    for(reference in references) {
      cat(paste0("Evaluating P=", p, ", PARTIAL=", partial, ", REF=", reference, "\n"))
      res <- dbGetQuery(conn, paste0("SELECT ",
                                     "datasets.UUID AS UUID, ",
                                     "METHOD, ",
                                     "PARTIAL_INFO, ",
                                     "BASELINE_TYPE, ",
                                     "datasets.BASELINE_CALLS AS ORACLE_BASELINE, ",
                                     "CALLS, ",
                                     "results.BASELINE_CALLS AS SELF_BASELINE, ",
                                     "P, ",
                                     "CORRP, ",
                                     "LOG_MEAN, ",
                                     "PERTURBATION, ",
                                     "REP_NOISE, ",
                                     "FC_ABSOLUTE, ",
                                     "FC_RELATIVE, ",
                                     "FC_PARTIAL, ",
                                     "MED_ABS_TOTAL, ",
                                     "MED_REL_TOTAL, ",
                                     "PERCENT_DIFF, ",
                                     "TPR, ",
                                     "FPR ",
                                     "FROM results LEFT JOIN datasets ON ",
                                     "results.UUID=datasets.UUID ",
                                     "WHERE P=", p, " ",
                                     "AND PARTIAL_INFO=", partial, " ",
                                     "AND BASELINE_TYPE='",reference,"';"))
      
      # Strip "result-less" entries
      res <- res %>%
        filter(!is.na(TPR) & !is.na(FPR))

      # Remove some of the simulations with huge median abundances
      # res <- res %>%
      #   filter(MED_ABS_TOTAL < 2e06)
      
      # Remove simulations where > 50% of features are differential
      # res <- res %>%
      #   filter(PERCENT_DIFF < 0.5)
      
      # Generate a fold change low/med/high factor labeling
      res$FC_plot <- sapply(res$FC_ABSOLUTE, function(x) {
        if(x < 1) {
          1 / x
        } else {
          x
        }
      })
      res$FC_plot <- cut(res$FC_plot, breaks = qq)
      levels(res$FC_plot) <- c("low", "moderate", "high")
      
      # Render median absolute totals on log scale for better visualization
      res$LOG_MAT <- log(res$MED_ABS_TOTAL)
      
      # Shuffle before plotting
      res <- res[sample(1:nrow(res)),]
      
      plot_tag <- paste0("p-", p, "_partial-", partial, "_ref-", reference)
      plot_ROC(res, plot_tag)
      plot_ROC(res %>% filter(PERCENT_DIFF <= 0.5), plot_tag = NULL, "FC_plot", "Fold change") # Look at no. features DE?
      plot_ROC(res, plot_tag = NULL, "LOG_MAT", "Log med.\ntotal (orig.)")
      # plot_ROC(res, plot_tag, "MED_REL_TOTAL", "Median total relative abundance")
      plot_ROC(res, plot_tag = NULL, "PERCENT_DIFF", "Percent DA features")
      # plot_ROC(res, plot_tag = NULL, "CORRP", "Correlation level")
      # plot_ROC(res, plot_tag = NULL, "LOG_MEAN", "Log mean abundance")
      # plot_ROC(res, plot_tag = NULL, "PERTURBATION", "Log mean perturbation")
      # plot_ROC(res, plot_tag, "REP_NOISE", "Noise within condition")
      
      # Does high abundance and high replicate noise yield big FPR?
      # label_combo(res, plot_tag, res$LOG_MAT > 17 & res$FC_plot == "high")
      # label_combo(res, plot_tag, res$LOG_MAT > 16 & res$REP_NOISE < 0.25)
      # label_combo(res, plot_tag, res$PERCENT_DIFF > 0.5)
      label_combo(res, plot_tag, res$MED_REL_TOTAL < 500000)
      
      # LM test - What explains the high FPR in some samples?
      # Answer - Surprisingly, it's having a high starting abundance. These are
      #          generally scenarios where you've got a big downsampling
      
      # lm_res <- res %>%
      #   select(TPR, FPR, METHOD, CORRP, LOG_MEAN, PERTURBATION, REP_NOISE, FC_ABSOLUTE, FC_RELATIVE, PERCENT_DIFF, MED_ABS_TOTAL) %>%
      #   mutate(FC_ABSOLUTE = log(FC_ABSOLUTE)) %>%
      #   mutate(MED_ABS_TOTAL = log(MED_ABS_TOTAL))
      # fit_fpr <- randomForest(FPR ~ ., data = lm_res %>% select(-TPR))
      # varImpPlot(fit_fpr)
      # ggplot(data.frame(x = lm_res$FPR, y = fit_fpr$predicted),
      #        aes(x = x, y = y, fill = lm_res$MED_ABS_TOTAL)) +
      #   geom_point(size = 2, shape = 21) +
      #   scale_fill_distiller(palette = "RdYlBu") + # continuous
      #   # scale_fill_brewer(palette = "RdYlBu") + # discrete
      #   labs(x = "True FPR",
      #        y = "Predicted FPR",
      #        fill = "MED_ABS_TOTAL")
      
      # Find "interesting" simulations and visualize them as stacked bar plots
      # temp <- res %>%
      #     filter(LOG_MAT > 16) %>%
      #     filter(FPR > 0.65)
      # 
      # temp2 <- temp[sample(1:nrow(temp), size = 1),]
      # str(temp2)
      # 
      # data <- readRDS(file.path("output", "datasets", paste0(temp2$UUID, ".rds")))
      # palette <- generate_highcontrast_palette(5000)
      # plot_stacked_bars(data$simulation$abundances, palette = palette)
      # plot_stacked_bars(data$simulation$observed_counts1, palette = palette)
      
      # What do false positive calls look like here vs. true negatives?
      # Like compositional effects.
      # But why are they more prevalent where original total abundances are 
      # high?
      # rates <- calc_DA_discrepancy(temp2$CALLS, temp2$ORACLE_BASELINE)
      # idx <- sample(which(rates$FP_calls == TRUE), size = 1)
      # par(mfrow = c(1, 2))
      # plot(data$simulation$abundances[,idx])
      # plot(data$simulation$observed_counts1[,idx])
    }
  }
}

# ------------------------------------------------------------------------------
#   Distributions of simulated datasets
# ------------------------------------------------------------------------------

res <- dbGetQuery(conn, "SELECT P, CORRP, PERCENT_DIFF, FC_ABSOLUTE FROM datasets")
res <- res %>%
  filter(CORRP %in% c(0, 4))
res$CORRP <- factor(res$CORRP, levels = c("0", "4"))
levels(res$CORRP) <- c("Independent features", "Strongly correlated features")
pl <- ggplot(res, aes(x = PERCENT_DIFF, fill = factor(P))) +
  geom_density(alpha = 0.6) +
  scale_fill_brewer(palette = "RdYlGn") +
  facet_wrap(. ~ CORRP) +
  xlim(c(0, 1)) +
  theme_bw() +
  labs(x = "percent differentially abundant features",
       fill = "Feature number")
show(pl)
ggsave(file.path("output",
                 "images",
                 paste0("PERCENT_DIFF.png")),
       plot = pl,
       dpi = 100,
       units = "in",
       height = 3,
       width = 6)

pl <- ggplot(res, aes(x = FC_ABSOLUTE, fill = factor(P))) +
  geom_density(alpha = 0.6) +
  scale_fill_brewer(palette = "RdYlGn") +
  facet_wrap(. ~ CORRP) +
  xlim(c(0, 10)) +
  theme_bw() +
  labs(x = "absolute fold change in abundance",
       fill = "Feature number")
show(pl)
ggsave(file.path("output",
                 "images",
                 paste0("FC_ABSOLUTE.png")),
       plot = pl,
       dpi = 100,
       units = "in",
       height = 3,
       width = 6)

# ------------------------------------------------------------------------------
#   Information "restored" in partial case
#
#   Whoops, it looks like the "partially informative" totals scenario 
#   essentially fully restores total abundances!
# ------------------------------------------------------------------------------

res <- dbGetQuery(conn, paste0("SELECT datasets.UUID, FC_ABSOLUTE, TOTALS_C_FC, PARTIAL FROM ",
                               "datasets LEFT JOIN characteristics ON ",
                               "datasets.UUID = characteristics.UUID WHERE ",
                               "PARTIAL = 1"))

# Make this symmetric, so it's comparable to TOTAL_C_FC
res$FC_ABSOLUTE <- sapply(res$FC_ABSOLUTE, function(x) {
  ifelse(x < 1, 1/x, x)
})

pl <- ggplot(res, aes(x = FC_ABSOLUTE, y = TOTALS_C_FC)) +
  geom_point(size = 2, shape = 21, fill = "#bbbbbb") +
  theme_bw() +
  labs(x = "Absolute fold change in totals (original abundances)",
       y = "Absolute fold change in totals (observed abundances)")
show(pl)
ggsave(file.path("output",
                 "images",
                 paste0("PARTIAL_correlation.png")),
       plot = pl,
       dpi = 100,
       units = "in",
       height = 5,
       width = 5.5)

cor(res$FC_ABSOLUTE, res$TOTALS_C_FC)**2

# ------------------------------------------------------------------------------
#   Other / TBD
# ------------------------------------------------------------------------------

dbDisconnect(conn)
