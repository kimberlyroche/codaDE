library(codaDE)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

getPalette <- colorRampPalette(c("#58D68D", "#FFFFFF", "#9B59B6"))
gv2 <- getPalette(2)
gv4 <- getPalette(4)

# ------------------------------------------------------------------------------
#   Parse and plot Klein et al. library size vs. Gapdh expression
#
#   This code has been copied from the `parse_Klein` function
# ------------------------------------------------------------------------------

counts_A <- read.table(file.path("data",
                                 "Klein_2015",
                                 "GSM1599494_ES_d0_main.csv",
                                 "GSM1599494_ES_d0_main.csv"),
                       sep = ",")
counts_B <- read.table(file.path("data",
                                 "Klein_2015",
                                 "GSM1599497_ES_d2_LIFminus.csv",
                                 "GSM1599497_ES_d2_LIFminus.csv"),
                       sep = ",")
# Gene IDs same order
gene_IDs <- counts_A$V1
counts_A <- counts_A[,2:ncol(counts_A)]
counts_B <- counts_B[,2:ncol(counts_B)]

gapdh_idx <- which(gene_IDs == "Gapdh")
counts <- cbind(counts_A, counts_B)
groups <- c(rep("unstimulated", ncol(counts_A)), rep("LIF-2hr", ncol(counts_B)))

# Remove samples with low counts
gapdh_expr <- unlist(counts[gapdh_idx,])
remove_idx <- which(colSums(counts) < 5000 | gapdh_expr < 1)
counts <- counts[,-remove_idx]
groups <- groups[-remove_idx]

gapdh_expr <- unlist(counts[gapdh_idx,])

x <- unname(unlist(colSums(counts)))
y <- gapdh_expr

plot_df_Klein <- data.frame(x = gapdh_expr,
                      y = unname(unlist(colSums(counts))),
                      group = groups)

p1 <- ggplot(plot_df_Klein, aes(x = x, y = y)) +
  geom_smooth(formula = "y ~ x", method = "lm", aes(group = group, color = group), alpha = 0.5) +
  geom_point(size = 2, shape = 21, mapping = aes(fill = group)) +
  scale_fill_manual(values = gv2) +
  scale_color_manual(values = gv2) +
  labs(x = "Gapdh abundance",
       y = "total abundance",
       fill = "Treatment") +
  theme_bw() +
  theme(axis.title.x = element_text(margin = ggplot2::margin(t = 0.15, unit = "in")),
        axis.title.y = element_text(margin = ggplot2::margin(r = 0.15, unit = "in")),
        plot.margin = unit(c(0, 0.3, 0, 0.05), "in"),
        legend.position = "bottom") +
  guides(color = FALSE)

# ------------------------------------------------------------------------------
#   Parse and plot Hashimshony et al. library size vs. spike-in expression
#
#   This code has been copied from the `parse_Hashimshony` function
# ------------------------------------------------------------------------------

# Parse data and assignments
data <- read.table(file.path("data",
                             "Hashimshony_2016",
                             "GSE78779_Expression_C1_96_cells.txt"),
                   sep = "\t",
                   header = TRUE)

groups <- sapply(colnames(data)[2:ncol(data)], function(x) {
  strsplit(x, "_")[[1]][2]
})
# Eliminate samples with missing/ambiguous cell cycle markers
retain_idx <- which(!(groups %in% c(".", "n.a.")))
data <- data[,c(1, retain_idx+1)]
groups <- groups[retain_idx]
groups <- factor(unname(groups), levels = c("0", "0.3", "0.5", "1"))

spike_idx <- which(sapply(data$Sample, function(x) str_detect(x, "^ERCC-")))
data <- data[,2:ncol(data)]

spike_expr <- data[spike_idx,]
counts <- data[-spike_idx,]

plot_df_Hashimshony <- data.frame(x = unname(unlist(apply(spike_expr, 2, mean))),
                                  y = unname(unlist(colSums(counts))),
                                  group = groups)
p2 <- ggplot(plot_df_Hashimshony, aes(x = x, y = y)) +
  geom_smooth(formula = "y ~ x", method = "lm", color = "#888888", alpha = 0.33) +
  geom_point(size = 2, shape = 21, mapping = aes(fill = group)) +
  scale_fill_manual(values = gv4) +
  labs(x = "mean spike-in abundance",
       y = "total abundance",
       fill = "cell cycle progress") +
  theme_bw() +
  theme(axis.title.x = element_text(margin = ggplot2::margin(t = 0.15, unit = "in")),
        axis.title.y = element_text(margin = ggplot2::margin(r = 0.15, unit = "in")),
        plot.margin = unit(c(0, 0, 0, 0.05), "in"),
        legend.position = "bottom")

# ------------------------------------------------------------------------------
#   Combine plots
# ------------------------------------------------------------------------------

p <- plot_grid(p1, p2, ncol = 2, labels = c("a", "b"),
               rel_widths = c(1, 0.9), scale = 0.95, label_size = 18)
show(p)
ggsave(file.path("output", "images", "library_sizes.png"),
       plot = p,
       dpi = 100,
       units = "in",
       height = 4.5,
       width = 9)
