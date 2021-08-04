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

# ------------------------------------------------------------------------------
#   Distributions of simulated datasets
# ------------------------------------------------------------------------------

res <- dbGetQuery(conn, paste0("SELECT P, CORRP, PERCENT_DIFF_REALIZ, ",
                               "FC_ABSOLUTE, FC_PARTIAL FROM datasets ",
                               "WHERE FC_ABSOLUTE <= 10 AND FC_ABSOLUTE >= 0.1"))

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
show(pl)
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
show(pl)
ggsave(file.path("output",
                 "images",
                 paste0("FC_ABSOLUTE.png")),
       plot = pl,
       dpi = 100,
       units = "in",
       height = 3,
       width = 6)

# pl <- ggplot(temp, aes(x = FC_PARTIAL, fill = factor(P))) +
#   geom_density(alpha = 0.6) +
#   scale_fill_brewer(palette = "RdYlGn") +
#   facet_wrap(. ~ CORRP) +
#   xlim(c(0, 10)) +
#   theme_bw() +
#   labs(x = "absolute fold change in abundance",
#        fill = "Feature number")
# show(pl)
# ggsave(file.path("output",
#                  "images",
#                  paste0("FC_PARTIAL.png")),
#        plot = pl,
#        dpi = 100,
#        units = "in",
#        height = 3,
#        width = 6)

# ------------------------------------------------------------------------------
#   Information "restored" in partial case
# ------------------------------------------------------------------------------

plot_df <- data.frame(x = res$FC_ABSOLUTE, y = res$FC_PARTIAL)

pl <- ggplot(plot_df, aes(x = x, y = y)) +
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

cor(plot_df$x, plot_df$y)**2

# ------------------------------------------------------------------------------
#   FPR as a function of "remaining" missing information
#
#   This is interesting. I think the size of the discrepancy being rescued by
#   these "partially informative" abundances is still just way too small to
#   help much with accuracy.
# ------------------------------------------------------------------------------

DE_method <- "scran"

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

res <- dbGetQuery(conn, paste0("SELECT PARTIAL_INFO, FC_ABSOLUTE, ",
                               "FC_PARTIAL, FPR ",
                               "FROM results LEFT JOIN datasets ",
                               "ON results.UUID = datasets.UUID ",
                               "WHERE FC_ABSOLUTE <= 10 ",
                               "AND FC_ABSOLUTE >= 0.1 ",
                               "AND METHOD='", DE_method, "' ",
                               "AND P=1000"))

dbDisconnect(conn)

# Subset to cases where absolute and partially accurate abundances register
# an increase in condition B *AND* partials don't accidentally create a larger
# discrepancy
temp <- res %>%
  filter(FC_ABSOLUTE > 1 & FC_ABSOLUTE > FC_PARTIAL)

# Calculate the absolute difference in fold change; too stupid?
temp$delta <- temp$FC_ABSOLUTE - temp$FC_PARTIAL
# temp$delta <- temp$FC_PARTIAL / temp$FC_ABSOLUTE

x <- c()
y <- c()
i <- 1
while(i < nrow(temp)) {
  fpr_before <- temp[i,]$FPR
  fpr_after <- temp[i+1,]$FPR
  x <- c(x, temp[i,]$delta)
  y <- c(y, fpr_after - fpr_before)
  i <- i + 2
}

pl <- ggplot(data.frame(x = x, y = y), aes(x = x, y = y)) +
  geom_point(size = 2, shape = 21, fill = "#bbbbbb") +
  theme_bw() +
  labs(x = "Difference in fold change: absolute - partially informative abundances",
       y = paste0("Change in FPR w/ partial info (", DE_method, ")"))
show(pl)
ggsave(file.path("output",
                 "images",
                 paste0("partial_rescue_", DE_method, ".png")),
       plot = pl,
       dpi = 100,
       units = "in",
       height = 5,
       width = 5.5)
