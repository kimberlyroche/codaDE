setwd("C:/Users/kimbe/Documents/codaDE")

source("path_fix.R")

library(tidyverse)
library(codaDE)
library(RSQLite)
library(gridExtra)
library(RColorBrewer)

source("ggplot_fix.R")

plot_ROC <- function(data, fill_var = NULL, fill_var_label = NULL) {
  if(is.null(fill_var)) {
    fill_var <- "METHOD"
    fill_var_label <- "Method"
  }
  data$fpr <- 1 - data$fpr
  if(fill_var == "flag") {
    # Layer the points
    p <- ggplot() +
      geom_point(data = data[data$flag == FALSE,],
                 mapping = aes_string(x = "fpr", y = "tpr"), size = 2, shape = 21, fill = "#dddddd") +
      geom_point(data = data[data$flag == TRUE,],
                 mapping = aes_string(x = "fpr", y = "tpr"), size = 3, shape = 21, fill = "#1ab079") +
      xlim(c(0,1)) +
      ylim(c(0,1)) +
      labs(x = "specificity (1 - FPR)",
           y = "sensitivity (TPR)",
           fill = fill_var_label) +
      theme_bw() +
      facet_wrap(. ~ METHOD, ncol = 4)
  } else {
    p <- ggplot(data, aes_string(x = "fpr", y = "tpr", fill = fill_var)) +
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
    p <- p +
      scale_fill_distiller(palette = "RdYlBu")
  }
  if(is.factor(data[[fill_var]])) {
    p <- p +
      scale_fill_brewer(palette = "RdYlBu")
  }
  if(fill_var == "METHOD") {
    p <- p +
      scale_fill_manual(values = palette)
  }
  # ggsave(file.path("output", "images", paste0(fill_var, ".png")),
  #        dpi = 100,
  #        units = "in",
  #        height = 3,
  #        width = 10)
  p
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
dir.create(file.path("output", "images", "full_results"), showWarnings = FALSE)
dir.create(file.path("output", "images", "full_results", "self_ref"), showWarnings = FALSE)
dir.create(file.path("output", "images", "full_results", "oracle_ref"), showWarnings = FALSE)
dir.create(file.path("output", "images", "full_results", "self_ref", "partial"), showWarnings = FALSE)
dir.create(file.path("output", "images", "full_results", "self_ref", "no_partial"), showWarnings = FALSE)
dir.create(file.path("output", "images", "full_results", "oracle_ref", "partial"), showWarnings = FALSE)
dir.create(file.path("output", "images", "full_results", "oracle_ref", "no_partial"), showWarnings = FALSE)

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

ps <- c(100, 1000, 5000)
partials <- c(0, 1)
references <- c("oracle", "self")
# methods <- c("ALDEx2", "DESeq2", "MAST", "scran")

# Build a fold change color scale
ramped_palette1 <- colorRampPalette(c("#00cc00", "#ffffff", "#fe8b00"))(7)
res <- dbGetQuery(conn, "SELECT FC_ABSOLUTE FROM datasets")$FC_ABSOLUTE
res <- sapply(res, function(x) {
  if(x < 1) {
    1 / x
  } else {
    x
  }
})
qq <- quantile(res, probs = seq(from = 0, to = 1, length.out = 4))
ramped_palette2 <- colorRampPalette(c("#7479c4", "#e63030"))(3)

for(p in ps) {
  for(partial in partials) {
    for(reference in references) {
      res <- dbGetQuery(conn, paste0("SELECT ",
                                     "datasets.UUID AS UUID, ",
                                     "METHOD, ",
                                     "PARTIAL_INFO, ",
                                     "BASELINE_TYPE, ",
                                     "datasets.BASELINE_CALLS AS ORACLE_BASELINE, ",
                                     "CALLS, ",
                                     "P, ",
                                     "CORRP, ",
                                     "LOG_MEAN, ",
                                     "PERTURBATION, ",
                                     "REP_NOISE, ",
                                     "FC_ABSOLUTE, ",
                                     "FC_RELATIVE, ",
                                     "FC_PARTIAL, ",
                                     "results.BASELINE_CALLS AS SELF_BASELINE ",
                                     "FROM results LEFT JOIN datasets ON ",
                                     "results.UUID=datasets.UUID WHERE ",
                                     "P=", p,
                                     # " AND CORRP=", corrp,
                                     " AND PARTIAL_INFO=", partial,
                                     # " AND METHOD='", method, "'",
                                     " AND BASELINE_TYPE='",reference,"'",
                                     " AND CALLS IS NOT NULL"))
      
      # Manually calculate TPR and FPR
      res$tpr <- NA
      res$fpr <- NA
      res$percent_diff <- NA
      for(i in 1:nrow(res)) {
        if(res[i,]$BASELINE_TYPE == "oracle") {
          rates <- calc_DA_discrepancy(res[i,]$CALLS, res[i,]$ORACLE_BASELINE)
        } else {
          rates <- calc_DA_discrepancy(res[i,]$CALLS, res[i,]$SELF_BASELINE)
        }
        res[i,]$tpr <- rates$TPR
        res[i,]$fpr <- rates$FPR
        res[i,]$percent_diff <- rates$percent_DA
      }
      str(res)
      
      # Generate a fold change labeling
      res$FC_plot <- sapply(res$FC_ABSOLUTE, function(x) {
        if(x < 1) {
          1 / x
        } else {
          x
        }
      })
      res$FC_plot <- cut(res$FC_plot, breaks = qq)
      levels(res$FC_plot) <- c("low", "moderate", "high")
      
      # Hard to know what role changes in sequencing depth plan in this...
      
      plot_ROC(res)
      plot_ROC(res, "FC_plot", "Fold change") # Look at no. features DE?
      plot_ROC(res, "percent_diff", "Percent DA features")
      plot_ROC(res, "CORRP", "Correlation level")
      plot_ROC(res, "LOG_MEAN", "Log mean abundance")
      plot_ROC(res, "PERTURBATION", "Log mean perturbation")
      plot_ROC(res, "REP_NOISE", "Noise within condition")
      
      # Find simulations where DESeq2 had huge FPR
      temp <- res %>%
        filter(METHOD == "DESeq2") %>%
        filter(!is.na(fpr)) %>%
        filter(fpr < 0.2) %>%
        filter(percent_diff > 0.5)
      str(temp)
      hist(temp$LOG_MEAN)
      
      temp2 <- temp[sample(1:nrow(temp), size = 1),]
      str(temp2)
      data <- readRDS(file.path("output", "datasets", paste0(temp2$UUID, ".rds")))
      
      res$flag <- FALSE
      # res$flag[res$FC_ABSOLUTE > 4 & res$LOG_MEAN >= 4 & res$CORRP >= 2] <- TRUE
      # res$flag[res$FC_plot == "low"] <- TRUE
      res$flag[res$UUID == temp2$UUID] <- TRUE
      plot_ROC(res, "flag", "selected")
      
      m1 <- mean(rowSums(data$simulation$abundances[1:10,]))
      m2 <- mean(rowSums(data$simulation$abundances[11:20,]))
      
      rates <- calc_DA_discrepancy(temp2$CALLS, temp2$ORACLE_BASELINE)
      fp_idx <- sample(which(rates$FP_calls == TRUE), size = 1)
      plot_df <- data.frame(x = 1:20,
                            y = data$simulation$abundances[,fp_idx],
                            condition = c(rep("A", 10), rep("B", 10)),
                            type = "abundances")
      plot_df <- rbind(plot_df,
                       data.frame(x = 1:20,
                                  y = data$simulation$observed_counts1[,fp_idx],
                                  condition = c(rep("A", 10), rep("B", 10)),
                                  type = "relative abundances"))
      ggplot(plot_df, aes(x = x, y = y, fill = factor(condition))) +
        geom_point(size = 2, shape = 21) +
        facet_wrap(. ~ type, ncol = 2) +
        labs(x = "sample index", y = "abundance", fill = "Condition")
      
      palette <- generate_highcontrast_palette(1000)
      plot_stacked_bars(data$simulation$abundances, palette = palette)
      plot_stacked_bars(data$simulation$observed_counts1, palette = palette)
      
      # ggsave(file.path("output",
      #                  "images",
      #                  "full_results",
      #                  ifelse(reference == "self",
      #                         "self_ref",
      #                         "threshold_ref"),
      #                  ifelse(partial == 0,
      #                         "no_partial",
      #                         "partial"),
      #                  paste0("SS_P",p,"_CORRP",corrp,".png")),
      #        plot = pl,
      #        dpi = 100,
      #        units = "in",
      #        height = 3,
      #        width = 12)
    }
  }
}

# ------------------------------------------------------------------------------
#   Distributions of simulated datasets
# ------------------------------------------------------------------------------

res <- dbGetQuery(conn, "SELECT P, CORRP, PERCENT_DIFF, FC_ABSOLUTE FROM datasets")
res$CORRP <- factor(res$CORRP)
levels(res$CORRP) <- c("independent feature simulations", "correlated feature simulations")

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
  xlim(c(0, 25)) +
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
#   Other / TBD
# ------------------------------------------------------------------------------

dbDisconnect(conn)









