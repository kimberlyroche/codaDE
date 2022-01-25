source("path_fix.R")

library(tidyverse)
library(codaDE)
library(cowplot)
library(RColorBrewer)
library(GGally)

# Create directories manually
dir.create("output", showWarnings = FALSE)
dir.create(file.path("output", "images"), showWarnings = FALSE)

# ------------------------------------------------------------------------------
#   How do real data features stack up with simulation for key characteristics?
#   Are poorly predicted datasets less typical with respect to informative
#     characteristics?
# ------------------------------------------------------------------------------

# result_files <- list.files(path = file.path("output",
#                                             "predictive_fits",
#                                             "oracle",
#                                             "regression_cpm",
#                                             "validation_results",
#                                             "no_norm"),
#                            pattern = "results_(.*?)_threshold1\\.tsv",
#                            full.names = TRUE)
# results <- NULL
# for(file in result_files) {
#   temp <- read.table(file, sep = "\t", header = TRUE)
#   results <- rbind(results,
#                    temp)
# }
# results <- results %>%
#   filter(result_type == "true") %>%
#   select(dataset, DE_method, score_type, point)
# results <- results %>%
#   pivot_wider(names_from = "score_type", values_from = "point")
# results$FPR <- 1 - results$FPR
# 
# datasets <- c("VieiraSilva", "Barlow", "Song", "Monaco", "Hagai", "Owens", "Klein", "Yu", "Muraro")
# thresholds <- rep(1, 9)
# 
# realdata_features <- NULL
# for(i in 1:length(datasets)) {
#   this_dataset <- datasets[i]
# 
#   cat(paste0("Parsing dataset ", this_dataset, "\n"))
#   # Parse real data
#   saved_data <- file.path("output",
#                           paste0("filtered_data_",this_dataset,"_threshold",thresholds[i],".rds"))
#   if(!file.exists(saved_data)) {
#     stop("Parsed data for ", this_dataset, " not found!\n")
#   }
#   data <- readRDS(saved_data)
# 
#   # This takes 2-3 min. to run on 15K features
#   features <- characterize_dataset(data$relative[data$groups == levels(data$groups)[1],],
#                                    data$relative[data$groups == levels(data$groups)[2],])
#   rm(data)
# 
#   realdata_features <- rbind(realdata_features,
#                              data.frame(dataset = this_dataset,
#                                         feature = names(features),
#                                         value = unname(unlist(features))))
# }
# 
# ggplot(realdata_features, aes(x = feature, y = value, fill = dataset)) +
#   geom_point(size = 2, shape = 21) +
#   theme_bw() +
#   coord_flip()

# Peek at the informative features
# k <- 3
# ta <- readRDS("output/predictive_fits/oracle/regression_cpm/oracle_TPR_ALDEx2_variable-importance.rds") %>%
#   arrange(desc(Overall)) %>%
#   slice(1:k)
# td <- readRDS("output/predictive_fits/oracle/regression_cpm/oracle_TPR_DESeq2_variable-importance.rds") %>%
#   arrange(desc(Overall)) %>%
#   slice(1:k)
# tx <- readRDS("output/predictive_fits/oracle/regression_cpm/oracle_TPR_scran_variable-importance.rds") %>%
#   arrange(desc(Overall)) %>%
#   slice(1:k)
# fa <- readRDS("output/predictive_fits/oracle/regression_cpm/oracle_FPR_ALDEx2_variable-importance.rds") %>%
#   arrange(desc(Overall)) %>%
#   slice(1:k)
# fd <- readRDS("output/predictive_fits/oracle/regression_cpm/oracle_FPR_DESeq2_variable-importance.rds") %>%
#   arrange(desc(Overall)) %>%
#   slice(1:k)
# fs <- readRDS("output/predictive_fits/oracle/regression_cpm/oracle_FPR_scran_variable-importance.rds") %>%
#   arrange(desc(Overall)) %>%
#   slice(1:k)
# 
# unique(c(rownames(ta), rownames(td), rownames(ts)))
# unique(c(rownames(fa), rownames(fd), rownames(fs)))

features <- readRDS(file.path("output", "feature_dump_cpm.rds"))

ffeat <- features %>%
  select(c(FW_LOG_SD_D, FW_CLR_PFC05_D, FW_CLR_PFC2_D, FPR))
ffeat$FPR <- 1 - ffeat$FPR

# # ------------------------------------------------------------------------------
# #   Render correlation plot
# # ------------------------------------------------------------------------------
# 
# # feature_name <- "FW_LOG_SD_D"
# # feature_name <- "FW_CLR_PFC05_D"
# feature_name <- "FW_CLR_PFC2_D"
# plot_piece1 <- ffeat %>%
#   select(c(feature_name, "FPR")) %>%
#   mutate(dataset = "simulation", DE_method = NA)
# plot_piece2 <- realdata_features %>%
#   filter(feature == feature_name) %>%
#   mutate(placeholder = value) %>%
#   select(-c(feature, value)) %>%
#   left_join(results, by = "dataset")
# colnames(plot_piece2)[2] <- feature_name
# 
# plot_df <- rbind(plot_piece1 %>% select(c("dataset", "DE_method", feature_name, "FPR")),
#                  plot_piece2 %>% select(c("dataset", "DE_method", feature_name, "FPR")))
# ggplot() +
#   geom_point(data = plot_piece1 %>% select(c("dataset", "DE_method", feature_name, "FPR")),
#              mapping = aes_string(x = feature_name, y = "FPR")) +
#   geom_point(data = plot_piece2 %>% select(c("dataset", "DE_method", feature_name, "FPR")),
#              mapping = aes_string(x = feature_name, y = "FPR", fill = "dataset"),
#              size = 3, shape = 21)

# ------------------------------------------------------------------------------
#   Render correlograms: FPR
# ------------------------------------------------------------------------------

colnames(ffeat) <- c("log count\nstandard deviation",
                     "percent strongly\ndecreasing CLR counts",
                     "1 - percent strongly\nincreasing CLR counts",
                     "specificity")

p <- ggpairs(ffeat, upper = list(continuous = wrap("cor", method = "spearman"))) +
  theme_bw()
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
