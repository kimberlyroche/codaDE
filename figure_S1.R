source("path_fix.R")

library(tidyverse)
library(codaDE)
library(RSQLite)
library(gridExtra)
library(RColorBrewer)
library(cowplot)
library(ggridges)

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

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

ramped_palette <- colorRampPalette(c("#7479c4", "#e63030"))(3)

res <- dbGetQuery(conn, paste0("SELECT P, CORRP, PERCENT_DIFF_REALIZ, ",
                               "FC_ABSOLUTE, FC_PARTIAL FROM datasets ",
                               "WHERE FC_ABSOLUTE <= 10 AND FC_ABSOLUTE >= 0.1"))

dbDisconnect(conn)

temp <- res %>%
  filter(CORRP %in% c(0, 4))
temp$CORRP <- factor(temp$CORRP, levels = c("0", "4"))
levels(temp$CORRP) <- c("Independent features", "Strongly correlated features")

lr_margin <- 0.45

# ------------------------------------------------------------------------------
#  Percent differential features
# ------------------------------------------------------------------------------

n_bins <- 30
bottom_offset <- 0.15

p1 <- ggplot(temp, aes(x = PERCENT_DIFF_REALIZ, y = factor(P), fill = factor(P))) +
  geom_density_ridges(stat = "binline", bins = n_bins, scale = 0.9) +
  scale_fill_brewer(palette = "RdYlGn") +
  scale_y_discrete(expand = expansion(add = c(bottom_offset, 1))) +
  xlim(c(0, 1)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = "percent differentially abundant features",
       fill = "Feature number") +
  theme(legend.position = "none",
        plot.margin = ggplot2::margin(t = 0.4, l = lr_margin, r = lr_margin, b = 0.5, "cm"))

# ------------------------------------------------------------------------------
#   Absolute fold change in totals
# ------------------------------------------------------------------------------

p2 <- ggplot(temp, aes(x = FC_ABSOLUTE, y = factor(P), fill = factor(P))) +
  geom_density_ridges(stat = "binline", bins = n_bins, scale = 0.9) +
  scale_fill_brewer(palette = "RdYlGn") +
  scale_y_discrete(expand = expansion(add = c(bottom_offset, 1))) +
  xlim(c(0, 10)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = "fold change in abundance",
       fill = "Feature number") +
  theme(legend.position = "none",
        plot.margin = ggplot2::margin(t = 0.4, l = lr_margin, r = lr_margin, b = 0.5, "cm"))

# ------------------------------------------------------------------------------
#   Percent zeros
# ------------------------------------------------------------------------------

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
res <- dbGetQuery(conn, paste0("SELECT P, COMP_C_P0_A, COMP_C_P0_B ",
                               "FROM characteristics LEFT JOIN datasets ",
                               "ON characteristics.UUID=datasets.UUID ",
                               "WHERE PARTIAL = 0;"))
dbDisconnect(conn)

temp2 <- data.frame(sparsity = c(res$COMP_C_P0_A, res$COMP_C_P0_B),
                    P = c(res$P, res$P))

p3 <- ggplot(temp2, aes(x = sparsity, y = factor(P), fill = factor(P))) +
  geom_density_ridges(stat = "binline", bins = n_bins, scale = 0.9) +
  scale_fill_brewer(palette = "RdYlGn") +
  scale_y_discrete(expand = expansion(add = c(bottom_offset, 1))) +
  xlim(c(0,1)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = "percent zeros",
       fill = "Feature number")
legend <- get_legend(p3 + theme(legend.position = "bottom"))
p3 <- p3  +
  theme(legend.position = "none",
        plot.margin = ggplot2::margin(t = 0.4, l = lr_margin, r = lr_margin, b = 0.5, "cm"))

prow <- plot_grid(p1,
                  p2,
                  p3,
                  align = 'vh',
                  hjust = -1,
                  nrow = 1)

pl <- plot_grid(prow,
                legend,
                ncol = 1,
                rel_heights = c(1, 0.1))

dev.off()
tiff(file.path("output", "images", "S1.tif"), units = "in", width = 10, height = 5, res = 300)
pl
dev.off()

# ggsave(file.path("output",
#                  "images",
#                  paste0("simulation_characteristics.svg")),
#        plot = pl,
#        dpi = 100,
#        units = "px",
#        height = 600,
#        width = 1000)
