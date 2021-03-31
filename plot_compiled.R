library(tidyverse)
library(codaDE)

setwd("C:/Users/kimbe/Documents/codaDE")

p <- 5000
rt <- "fpr"

data <- readRDS(file.path("output", paste0("simresults_p",p,"_simulated_all.rds")))

sub_data <- data %>%
  filter(rate_type == rt) %>%
  filter(delta_mean_v2 < 10)

if(!exists("cpalette")) {
  cpalette <- generate_highcontrast_palette(length(unique(sub_data$method)))
}

sub_data$method <- as.factor(sub_data$method)
# levels(sub_data$method) <- c("NB GLM (baseline)",
#                              "DESeq2",
#                              "edgeR (w/o TMM)",
#                              "MAST",
#                              "scran (findMarkers)",
#                              "simulated spike-in normalization",
#                              "Wilcox (via Seurat)")
levels(sub_data$method) <- c("NB GLM (baseline)",
                             "DESeq2")

pl <- ggplot(sub_data, aes(x = delta_mean_v2, y = rate, color = method)) +
  scale_color_manual(values = cpalette) +
  geom_point(alpha = 0.5) +
  geom_smooth() +
  labs(color = "method")
if(rt == "tpr") {
  pl <- pl +
    ylim(c(0.35, 1))
} else {
  pl <- pl +
    ylim(c(0, 0.3))
}
pl <- pl +
  xlab("change in mean abundance between conditions") +
  ylab("false positive rate")

show(pl)
ggsave(file.path("output", "images", paste0("DE_p",p,"_all_models_",rt,".png")),
       pl,
       dpi = 300,
       units = "in",
       height = 6.5,
       width = 9)
