library(tidyverse)
library(codaDE)

# ------------------------------------------------------------------------------
#   Figure 2 - Simulated absolute vs. observed abundances
# ------------------------------------------------------------------------------

# Stacked bar plots

sim_data <- simulate_sequence_counts(n = 5, p = 100, k = 1, ref_data = "simulated_bulk",
                                     asymmetry = 0.75, proportion_da = 0.6,
                                     spike_in = FALSE, possible_fold_changes = NULL)

palette <- generate_highcontrast_palette(ncol(sim_data$abundances))

plot_stacked_bars(sim_data$abundances, palette = palette, save_name = NULL)
ggsave(file.path("output", "images", "poster_F2a.png"),
       units = "in",
       dpi = 300,
       height = 6,
       width = 6)

plot_stacked_bars(sim_data$observed_counts1, palette = palette, save_name = NULL)
ggsave(file.path("output", "images", "poster_F2b.png"),
       units = "in",
       dpi = 300,
       height = 6,
       width = 6)

# ------------------------------------------------------------------------------
#   OMITTED: Athanasiadou vs. Barlow (or Morton) totals
# ------------------------------------------------------------------------------

ciona <- readRDS("data/absolute_Athanasiadou_ciona.rds")
palette <- generate_highcontrast_palette(100)
palette <- sample(palette, size = length(unique(ciona$groups)))

total_abundance <- as.data.frame(colSums(ciona$counts))
total_abundance <- cbind(rownames(total_abundance), ciona$groups, total_abundance)
colnames(total_abundance) <- c("sample", "cell_marker", "total_counts")
rownames(total_abundance) <- NULL

ggplot(total_abundance, aes(x = sample, y = total_counts, fill = cell_marker)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette) +
  theme_bw()
ggsave(file.path("output", "images", "poster_F2a.png"),
       units = "in",
       dpi = 300,
       height = 6,
       width = 4)

microbes <- readRDS("data/absolute_Morton.rds")
microbes <- readRDS("data/absolute_Barlow.rds")

palette <- generate_highcontrast_palette(100)
palette <- sample(palette, size = length(unique(microbes$groups)))

total_abundance <- as.data.frame(colSums(microbes$counts))
total_abundance <- cbind(rownames(total_abundance), microbes$groups, total_abundance)
colnames(total_abundance) <- c("sample", "event", "total_counts")
rownames(total_abundance) <- NULL
# Fix annoying ggplot x-axis alphabetizing
total_abundance$sample <- factor(total_abundance$sample, levels = total_abundance$sample)

ggplot(total_abundance, aes(x = sample, y = total_counts, fill = event)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette) +
  theme_bw()
ggsave(file.path("output", "images", "poster_F2b.png"),
       units = "in",
       dpi = 300,
       height = 6,
       width = 8)

# ------------------------------------------------------------------------------
#   FIGURES 4-5 - TBD (simulations running...)
# ------------------------------------------------------------------------------






