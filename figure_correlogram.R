source("path_fix.R")

library(tidyverse)
library(codaDE)
library(cowplot)
library(RColorBrewer)
library(GGally)

# Create directories manually
dir.create("output", showWarnings = FALSE)
dir.create(file.path("output", "images"), showWarnings = FALSE)

features <- readRDS(file.path("output", "feature_dump_cpm.rds"))

ffeat <- features %>%
  select(c(FW_LOG_SD_D, FW_CLR_PFC05_D, FW_CLR_PFC2_D, FPR))
ffeat$FPR <- 1 - ffeat$FPR

# ------------------------------------------------------------------------------
#   Render correlograms: FPR
# ------------------------------------------------------------------------------

colnames(ffeat) <- c("log count\nstandard deviation",
                     "percent strongly\ndecreasing CLR counts",
                     "1 - percent strongly\nincreasing CLR counts",
                     "specificity")

p <- ggpairs(ffeat, upper = list(continuous = wrap("cor", method = "spearman"))) +
  theme_bw()
p
ggsave("output/images/FPR_correl.png",
       p,
       dpi = 100,
       units = "in",
       height = 6.5,
       width = 6.5)

# ------------------------------------------------------------------------------
#   Render correlograms: TPR
# ------------------------------------------------------------------------------

tfeat <- features %>%
  select(c(CORR_RA_SKEW, COMP_C_P0_A, COMP_C_P0_B, TPR))
colnames(tfeat) <- c("skew in correlation\nof relative abundances",
                     "percent zeros\ncondition A",
                     "percent zeros\ncondition B",
                     "sensitivity")

p <- ggpairs(tfeat, upper = list(continuous = wrap("cor", method = "spearman"))) +
  theme_bw()
ggsave("output/images/TPR_correl.png",
       p,
       dpi = 100,
       units = "in",
       height = 6.5,
       width = 6.5)
