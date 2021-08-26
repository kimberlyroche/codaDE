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
                 size = 3,
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

qq_res <- dbGetQuery(conn, paste0("SELECT FC_ABSOLUTE FROM datasets ",
                                  "WHERE FC_ABSOLUTE <= 10 AND ",
                                  "FC_ABSOLUTE >= 0.1"))$FC_ABSOLUTE
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
                                     "PERCENT_DIFF_REALIZ, ",
                                     "TPR, ",
                                     "FPR ",
                                     "FROM results LEFT JOIN datasets ON ",
                                     "results.UUID=datasets.UUID ",
                                     "WHERE P=", p, " ",
                                     "AND PARTIAL_INFO=", partial, " ",
                                     "AND BASELINE_TYPE='",reference,"' ",
                                     "AND FC_ABSOLUTE <= 10 ",
                                     "AND FC_ABSOLUTE >= 0.1;"))
      
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
      
      # Shuffle before plotting
      res <- res[sample(1:nrow(res)),]
      
      # plot_tag <- paste0("p-", p, "_partial-", partial, "_ref-", reference)
      # plot_ROC(res, plot_tag, "FC_plot", "Fold change") # Look at no. features DE?

      res$sign <- ifelse(res$FC_ABSOLUTE > 1, 1, -1)
      plot_ROC(res %>%
                 filter(FC_plot == "high") %>%
                 filter(FPR < 0.05 | FPR > 0.5) %>%
                 filter(TPR > 0.8), plot_tag = NULL, "sign")
      plot_ROC(res %>%
                 filter(FC_plot == "high") %>%
                 filter(FPR < 0.05 | FPR > 0.5) %>%
                 filter(TPR > 0.8), plot_tag = NULL, "LOG_MEAN")
      
      plot_tag <- paste0("p-", p, "_partial-", partial, "_ref-", reference, "_max50")
      plot_ROC(res %>% filter(PERCENT_DIFF_REALIZ < 0.5), plot_tag, "FC_plot", "Fold change") # Look at no. features DE?
      
      # plot_ROC(res, plot_tag)
      # plot_ROC(res, plot_tag, "FC_plot", "Fold change") # Look at no. features DE?
      # plot_ROC(res, plot_tag = NULL, "LOG_MAT", "Log med.\ntotal (orig.)")
      # plot_ROC(res, plot_tag, "MED_REL_TOTAL", "Median total relative abundance")
      # plot_ROC(res, plot_tag = NULL, "PERCENT_DIFF", "Percent DA features")
      # plot_ROC(res, plot_tag = NULL, "CORRP", "Correlation level")
      # plot_ROC(res, plot_tag = NULL, "LOG_MEAN", "Log mean abundance")
      # plot_ROC(res, plot_tag = NULL, "PERTURBATION", "Log mean perturbation")
      # plot_ROC(res, plot_tag, "REP_NOISE", "Noise within condition")
      
      # Does high abundance and high replicate noise yield big FPR?
      # label_combo(res, plot_tag, res$LOG_MAT > 17 & res$FC_plot == "high")
      # label_combo(res, plot_tag, res$LOG_MAT > 16 & res$REP_NOISE < 0.25)
      # label_combo(res, plot_tag, res$PERCENT_DIFF > 0.5)
      # label_combo(res, plot_tag, res$PERCENT_DIFF < 0.5)
      
      # LM test - What explains the high FPR in some samples?
      # Answer - Surprisingly, it's having a high starting abundance. These are
      #          generally scenarios where you've got a big downsampling
      
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
      
      res <- res %>%
        filter(METHOD == 'DESeq2' & FPR > 0.6)
      uuid <- sample(res$UUID, size = 1) # Pull this dataset

      res2 <- dbGetQuery(conn, paste0("SELECT datasets.UUID, BASELINE_TYPE, CALLS, ",
                                      "datasets.BASELINE_CALLS AS ORACLE_CALLS, ",
                                      "results.BASELINE_CALLS AS SELF_CALLS ",
                                      "FROM results LEFT JOIN datasets ",
                                      "ON results.UUID=datasets.UUID WHERE ",
                                      "METHOD='DESeq2' AND ",
                                      "datasets.UUID='", uuid, "' AND ",
                                      "PARTIAL_INFO=0"))

      oracle_calls <- as.numeric(str_split(res2 %>%
                                             filter(is.na(SELF_CALLS)) %>%
                                             pull(ORACLE_CALLS), ";")[[1]])
      self_calls <- as.numeric(str_split(res2 %>%
                                           filter(!is.na(SELF_CALLS)) %>%
                                           pull(SELF_CALLS), ";")[[1]])

      oracle_adj <- p.adjust(oracle_calls, method = "BH")
      self_adj <- p.adjust(self_calls, method = "BH")

      idx <- which(oracle_adj >= 0.05 & self_adj < 0.05)

      data <- readRDS(file.path("output", "datasets", paste0(uuid, ".rds")))
      abundances <- data$simulation$abundances

      # What does a random disagreement (NO by oracle, YES by DESeq2) look like?
      ggplot(data.frame(x = 1:nrow(abundances),
                        y = abundances[,sample(idx, size = 1)],
                        group = factor(data$simulation$groups)),
             aes(x = x, y = y, fill = group)) +
        geom_point(size = 2, shape = 21) +
        theme_bw() +
        labs(x = "sample index",
             y = "abundance")
      
      # rates <- calc_DA_discrepancy(temp2$CALLS, temp2$ORACLE_BASELINE)
      # idx <- sample(which(rates$FP_calls == TRUE), size = 1)
      # par(mfrow = c(1, 2))
      # plot(data$simulation$abundances[,idx])
      # plot(data$simulation$observed_counts1[,idx])
    }
  }
}

dbDisconnect(conn)
