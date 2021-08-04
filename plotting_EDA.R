source("path_fix.R")

library(tidyverse)
library(codaDE)
library(RSQLite)
library(gridExtra)
library(RColorBrewer)

source("ggplot_fix.R")

# Create directories manually
dir.create("output", showWarnings = FALSE)
dir.create(file.path("output", "images"), showWarnings = FALSE)

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

qq_res <- dbGetQuery(conn, "SELECT FC_ABSOLUTE FROM datasets")$FC_ABSOLUTE
qq_res <- sapply(qq_res, function(x) {
  if(x < 1) {
    1 / x
  } else {
    x
  }
})
qq <- quantile(qq_res, probs = seq(from = 0, to = 1, length.out = 4))

ramped_palette <- colorRampPalette(c("#7479c4", "#e63030"))(3)

# ------------------------------------------------------------------------------
#   Distributions of simulated datasets
# ------------------------------------------------------------------------------

res <- dbGetQuery(conn, "SELECT P, CORRP, PERCENT_DIFF_REALIZ, FC_ABSOLUTE, FC_PARTIAL FROM datasets")

dbDisconnect(conn)

temp <- res %>%
  filter(CORRP %in% c(0, 4))
temp$CORRP <- factor(temp$CORRP, levels = c("0", "4"))
levels(temp$CORRP) <- c("Independent features", "Strongly correlated features")
pl <- ggplot(temp, aes(x = PERCENT_DIFF_REALIZ, fill = factor(P))) +
  geom_density(alpha = 0.6) +
  scale_fill_brewer(palette = "RdYlGn") +
  facet_wrap(. ~ CORRP) +
  xlim(c(0, 1)) +
  theme_bw() +
  labs(x = "percent differentially abundant features",
       fill = "Feature number")
# show(pl)
ggsave(file.path("output",
                 "images",
                 paste0("PERCENT_DIFF.png")),
       plot = pl,
       dpi = 100,
       units = "in",
       height = 3,
       width = 6)

pl <- ggplot(temp, aes(x = FC_ABSOLUTE, fill = factor(P))) +
  geom_density(alpha = 0.6) +
  scale_fill_brewer(palette = "RdYlGn") +
  facet_wrap(. ~ CORRP) +
  xlim(c(0, 10)) +
  theme_bw() +
  labs(x = "absolute fold change in abundance",
       fill = "Feature number")
# show(pl)
ggsave(file.path("output",
                 "images",
                 paste0("FC_ABSOLUTE.png")),
       plot = pl,
       dpi = 100,
       units = "in",
       height = 3,
       width = 6)

pl <- ggplot(temp, aes(x = FC_PARTIAL, fill = factor(P))) +
  geom_density(alpha = 0.6) +
  scale_fill_brewer(palette = "RdYlGn") +
  facet_wrap(. ~ CORRP) +
  xlim(c(0, 10)) +
  theme_bw() +
  labs(x = "absolute fold change in abundance",
       fill = "Feature number")
# show(pl)
ggsave(file.path("output",
                 "images",
                 paste0("FC_PARTIAL.png")),
       plot = pl,
       dpi = 100,
       units = "in",
       height = 3,
       width = 6)

# ------------------------------------------------------------------------------
#   Information "restored" in partial case
# ------------------------------------------------------------------------------

# res$FC_ABSOLUTE <- sapply(res$FC_ABSOLUTE, function(x) {
#   ifelse(x < 1, 1/x, x)
# })

pl <- ggplot(res, aes(x = FC_ABSOLUTE, y = FC_PARTIAL)) +
  geom_point(size = 2, shape = 21, fill = "#bbbbbb") +
  theme_bw() +
  labs(x = "Fold change in totals (original abundances)",
       y = "Fold change in totals (observed abundances)")
show(pl)
ggsave(file.path("output",
                 "images",
                 paste0("PARTIAL_CORR.png")),
       plot = pl,
       dpi = 100,
       units = "in",
       height = 5,
       width = 5.5)

cor(res$FC_ABSOLUTE, res$FC_PARTIAL)**2
