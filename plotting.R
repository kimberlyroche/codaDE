source("path_fix.R")

library(tidyverse)
library(codaDE)
library(RSQLite)
library(gridExtra)
library(RColorBrewer)
library(randomForest)

source("ggplot_fix.R")

string_to_calls <- function(call_string) {
  pvals <- as.numeric(strsplit(call_string, ";")[[1]])
  pvals <- p.adjust(pvals, method = "BH")
  pvals < 0.05
}

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
                 # size = 2,
                 # shape = 21,
                 color = "#dddddd") +
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
    pl <- ggplot() +
      # geom_point(data = data %>% filter(MED_ABS_TOTAL > 2e06),
      #            mapping = aes(x = FPR, y = TPR),
      #            color = "#dddddd") +
      # geom_point(data = data %>% filter(MED_ABS_TOTAL <= 2e06),
      #            mapping = aes_string(x = "FPR", y = "TPR", fill = fill_var),
      #            size = 2,
      #            shape = 21) +
      geom_point(data,
                 mapping = aes_string(x = "FPR", y = "TPR", fill = fill_var),
                 size = 2,
                 shape = 21) +
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
  show(pl)
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

plot_density <- function(data, plot_tag = NULL) {
  fill_var <- "METHOD"
  fill_var_label <- "Method"

  pl <- ggplot() +
    geom_density(data,
               mapping = aes(x = FPR, fill = METHOD)) +
    labs(x = "false discovery rate",
         y = "density",
         fill = fill_var_label) +
    scale_fill_manual(values = palette) +
    theme_bw() +
    theme(legend.position = "none") +
    facet_wrap(. ~ METHOD, ncol = 4, scales = "free_y")
  if(!is.null(plot_tag)) {
    ggsave(file.path("output", "images", paste0("FDR_",
                                                fill_var,
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
partials <- c(0)
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
qq <- c(qq[1:3], 5, Inf)

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
      
      # What proportion of calls exceed an FDR of 5% here?
      # res <- res %>%
      #   filter(METHOD == 'DESeq2') %>%
      #   mutate(FDR_big = ifelse(FPR > 0.05, "yes", "no")) %>%
      #   select(!c("ORACLE_BASELINE", "SELF_BASELINE", "CALLS"))
      # 
      # sum(res$FDR_big == "yes")/nrow(res)
      
      # Filter out what DESeq2 calls > 85% DE
      # res$DESeq2_DE <- sapply(res$SELF_BASELINE, function(x) sum(string_to_calls(x)))
      # res$DESeq2_DE <- res$DESeq2_DE / p
      # res <- res %>%
      #   filter(DESeq2_DE <= 0.85)
      
      # res <- res %>%
      #   filter(PERCENT_DIFF_REALIZ < 0.85)
      
      # ------------------------------------------------------------------------
      #   Histogram with real data overlaid
      # ------------------------------------------------------------------------
      
      # real_data <- readRDS("bad_good_results.rds")
      # real_data <- real_data %>%
      #   filter(dtype != "simulated") %>%
      #   select(dtype, fpr)
      # 
      # real_df <- data.frame(x = real_data$fpr,
      #                       y = rep(c(25,45), 4),
      #                       type = real_data$dtype)
      # ggplot() +
      #   geom_histogram(data = res[res$METHOD == "DESeq2",],
      #                  mapping = aes(x = FPR),
      #                  color = "white") +
      #   geom_point(data = real_df,
      #              mapping = aes(x = x, y = y, fill = type),
      #              size = 4,
      #              shape = 21,
      #              stroke = 1.2) +
      #   scale_fill_manual(values = generate_highcontrast_palette(8)) +
      #   labs(fill = "Data set") +
      #   theme_bw()
      
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
      levels(res$FC_plot) <- c("low", "moderate", "high", "very high")
      
      # Shuffle before plotting
      res <- res[sample(1:nrow(res)),]
      
      plot_tag <- paste0("p-", p, "_partial-", partial, "_ref-", reference)
      plot_ROC(res, plot_tag, "FC_plot", "Fold change")
      plot_density(res, plot_tag)
      
      # plot_ROC(res, plot_tag, "PERCENT_DIFF_REALIZ", "Realized proportion DE") # Look at no. features DE?
      
      # label_combo(res, plot_tag, res$PERCENT_DIFF_REALIZ < 0.5 & res$FC_plot %in% c("high","very high"))
    }
  }
}

dbDisconnect(conn)
