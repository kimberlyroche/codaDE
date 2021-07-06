setwd("C:/Users/kimbe/Documents/codaDE")

source("path_fix.R")

library(tidyverse)
library(codaDE)
library(RSQLite)
# library(ggROC)
# library(pROC)
library(gridExtra)
library(RColorBrewer)

source("ggplot_fix.R")

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
dir.create(file.path("output", "images", "full_results", "threshold_ref"), showWarnings = FALSE)
dir.create(file.path("output", "images", "full_results", "self_ref", "partial"), showWarnings = FALSE)
dir.create(file.path("output", "images", "full_results", "self_ref", "no_partial"), showWarnings = FALSE)
dir.create(file.path("output", "images", "full_results", "threshold_ref", "partial"), showWarnings = FALSE)
dir.create(file.path("output", "images", "full_results", "threshold_ref", "no_partial"), showWarnings = FALSE)

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

ps <- c(100, 1000, 5000)
corrps <- c(0, 1)
partials <- c(0, 1)
references <- c("threshold", "self")
# methods <- c("ALDEx2", "DESeq2", "MAST", "NBGLM", "scran")

# Build a fold change color scale
res <- dbGetQuery(conn, "SELECT FC_ABSOLUTE FROM datasets")$FC_ABSOLUTE
ramped_palette1 <- colorRampPalette(c("#00cc00", "#ffffff", "#fe8b00"))(7)
res <- dbGetQuery(conn, "SELECT FC_ABSOLUTE FROM datasets")$FC_ABSOLUTE
qq <- quantile(res, probs = seq(from = 0, to = 1, length.out = 4))
ramped_palette2 <- colorRampPalette(c("#7479c4", "#e63030"))(3)

for(p in ps) {
  for(corrp in corrps) {
    for(partial in partials) {
      for(reference in references) {
        res <- dbGetQuery(conn, paste0("SELECT * FROM results LEFT JOIN datasets ON ",
                                       "results.UUID=datasets.UUID WHERE ",
                                       "P=", p,
                                       " AND CORRP=", corrp,
                                       " AND PARTIAL_INFO=", partial,
                                       # " AND METHOD='", method, "'",
                                       " AND BASELINE='",reference,"'",
                                       " AND RESULT IS NOT NULL"))
        
        res <- res %>%
          pivot_wider(names_from = "RESULT_TYPE", values_from = "RESULT")
        
        # Plot with percent-differential labeling
        # PERCENT <- res$PERCENT_DIFF
        # PERCENT_factor <- cut(PERCENT, breaks = seq(from = 0.1, to = 0.8, by = 0.1))
        # # levels(FC_factor) <- c("low", "moderate", "high")
        # pl <- ggplot(res, aes(x = 1 - fpr, y = tpr, fill = PERCENT_factor)) +
        #   geom_point(size = 3, shape = 21) +
        #   scale_fill_manual(values = ramped_palette1) +
        #   xlim(c(0,1)) +
        #   ylim(c(0,1)) +
        #   labs(x = "specificity (1 - FPR)",
        #        y = "sensitivity (TPR)",
        #        fill = "Fold change ") +
        #   theme_bw() +
        #   theme(legend.position = "none") +
        #   facet_wrap(. ~ METHOD, ncol = 5)
        # show(pl)
        
        # Plot fold change labeling
        FC <- res$FC_ABSOLUTE
        FC_factor <- cut(FC, breaks = qq)
        levels(FC_factor) <- c("low (1-1.5x)", "moderate (1.5-2.5x)", "high (>2.5x)")
        pl <- ggplot(res, aes(x = 1 - fpr, y = tpr, fill = FC_factor)) +
          geom_point(size = 3, shape = 21) +
          scale_fill_manual(values = ramped_palette2) +
          xlim(c(0,1)) +
          ylim(c(0,1)) +
          labs(x = "specificity (1 - FPR)",
               y = "sensitivity (TPR)",
               fill = "Fold change ") +
          theme_bw() +
          # theme(legend.position = "none") +
          facet_wrap(. ~ METHOD, ncol = 5)
        # show(pl)
        ggsave(file.path("output",
                         "images",
                         "full_results",
                         ifelse(reference == "self",
                                "self_ref",
                                "threshold_ref"),
                         ifelse(partial == 0,
                                "no_partial",
                                "partial"),
                         paste0("SS_P",p,"_CORRP",corrp,"_labeled-FC.png")),
               plot = pl,
               dpi = 100,
               units = "in",
               height = 3,
               width = 14)
        
        # Plot with no labeling
        pl <- ggplot(res, aes(x = 1 - fpr, y = tpr, fill = METHOD)) +
          geom_point(size = 3, shape = 21) +
          scale_fill_manual(values = palette) +
          xlim(c(0,1)) +
          ylim(c(0,1)) +
          labs(x = "specificity (1 - FPR)",
               y = "sensitivity (TPR)",
               fill = "Method") +
          theme_bw() +
          theme(legend.position = "none") +
          facet_wrap(. ~ METHOD, ncol = 5)
        # show(pl)
        ggsave(file.path("output",
                         "images",
                         "full_results",
                         ifelse(reference == "self",
                                "self_ref",
                                "threshold_ref"),
                         ifelse(partial == 0,
                                "no_partial",
                                "partial"),
                         paste0("SS_P",p,"_CORRP",corrp,".png")),
               plot = pl,
               dpi = 100,
               units = "in",
               height = 3,
               width = 12)
      }
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









