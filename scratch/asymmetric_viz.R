source("path_fix.R")

library(tidyverse)
library(ggridges)
library(RColorBrewer)

# datasets <- list.files(file.path("output", "datasets"))
# results <- NULL
# for(i in 1:length(datasets)) {
#   if(i %% 100 == 0) {
#     cat(paste0("Iteration ", i, " / ", length(datasets), "\n"))
#   }
#   data <- readRDS(file.path("output", "datasets", datasets[i]))
#   a <- colMeans(data$simulation$abundances[1:10,])
#   b <- colMeans(data$simulation$abundances[11:20,])
#   results <- rbind(results,
#                    data.frame(file = datasets[i],
#                               P = length(a),
#                               percent_positive = sum(sign(b - a) > 0)/length(a)))
# }
# saveRDS(results, "asymmetry.rds")

results <- readRDS("asymmetry.rds")

ppalette <- colorRampPalette(brewer.pal(3, "Reds"))(3)

p <- ggplot(results, aes(x = percent_positive*100, y = factor(P), fill = factor(P))) +
  geom_density_ridges(stat = "binline", bins = 40, scale = 0.95) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = ppalette) +
  scale_y_discrete(expand = expansion(add = c(0.2, 1.0))) +
  labs(x = "\npercent features increasing across simulated conditions",
       y = "feature number\n")

ggsave(file.path("output", "images", "asymmetry.svg"),
       p,
       dpi = 100,
       units = "in",
       height = 5,
       width = 6)
