library(codaDE)
library(tidyverse)

palette <- generate_highcontrast_palette(8)
stats <- readRDS("output/real_data_filtering_statistics.rds")

plot_df <- stats %>%
  pivot_longer(!c(dataset, filter), names_to = "characteristic", values_to = "value")
plot_df$characteristic <- factor(plot_df$characteristic,
                                 levels = c("n_features",
                                            "stable_abs",
                                            "stable_rel",
                                            "zeros_abs",
                                            "zeros_rel"))
levels(plot_df$characteristic) <- c("Number of features",
                                    "Proportion of stable\nfeatures (absolute counts)",
                                    "Proportion of stable\nfeatures (relative abundances)",
                                    "Proportion of zeros\n(absolute counts)",
                                    "Proportion of zeros\n(relative abundances)")

ggplot(plot_df %>% filter(filter > 0),
       aes(x = filter, y = value, fill = dataset)) +
  geom_line() +
  geom_point(size = 2, shape = 21) +
  facet_wrap(. ~ dataset + characteristic, scales = "free_y", ncol = 5) +
  theme_bw() +
  labs(fill = "Dataset name") +
  scale_fill_manual(values = palette)
ggsave(file.path("output", "images", "realdata_characteristics.png"),
       dpi = 100,
       units = "in",
       height = 16,
       width = 14)
