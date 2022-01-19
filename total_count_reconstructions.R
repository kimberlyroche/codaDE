library(codaDE)
library(tidyverse)
library(cowplot)

datasets <- c("VieiraSilva", "Barlow", "Song", "Monaco", "Hagai", "Owens", "Klein", "Yu")
thresholds <- c(1, 1, 1, 1, 1, 1, 1, 1)

plots <- list()
for(i in 1:length(datasets)) {
  all_data <- readRDS(file.path("output", paste0("filtered_data_", datasets[i], "_threshold", thresholds[i], ".rds")))
  scaled_counts_A <- scaled_counts_ALDEx2(all_data$relative, pseudocount = 0.5)
  scaled_counts_D <- scaled_counts_DESeq2(all_data$relative, all_data$groups, pseudocount = 0.5)
  scaled_counts_S <- scaled_counts_scran(all_data$relative, all_data$groups, pseudocount = 0.5)
  
  plot(apply(all_data$relative, 1, function(x) mean(log(x + 0.5))))
  
  # plot_df <- data.frame(sample = 1:nrow(all_data$relative),
  #                       total = rowSums(all_data$relative),
  #                       group = all_data$groups,
  #                       type = "relative")
  plot_df <- data.frame(sample = 1:nrow(scaled_counts_A),
                        total = rowSums(scaled_counts_A),
                        group = all_data$groups,
                        type = "ALDEx2")
  plot_df <- rbind(plot_df,
                   data.frame(sample = 1:nrow(scaled_counts_D),
                              total = rowSums(scaled_counts_D),
                              group = all_data$groups,
                              type = "DESeq2"))
  plot_df <- rbind(plot_df,
                   data.frame(sample = 1:nrow(scaled_counts_S),
                              total = rowSums(scaled_counts_S),
                              group = all_data$groups,
                              type = "scran"))
  plot_df <- rbind(plot_df,
                   data.frame(sample = 1:nrow(all_data$absolute),
                              total = rowSums(all_data$absolute),
                              group = all_data$groups,
                              type = "absolute"))
  
  plot_df$type <- factor(plot_df$type, levels = c("relative",
                                                  "absolute",
                                                  "ALDEx2",
                                                  "DESeq2",
                                                  "scran"))
  
  plots[[i]] <- ggplot(plot_df, aes(x = sample, y = total, fill = group)) +
    geom_bar(stat = "identity") +
    facet_wrap(. ~ type, scales = "free_y", ncol = 4) +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2")
    # theme(legend.position = "none")
}

p <- plot_grid(plotlist = plots, ncol = 1)

ggsave("output/temp.png",
       p,
       dpi = 100,
       units = "in",
       height = 12,
       width = 8)

# Peak at false positives

data <- readRDS("output/filtered_data_Owens_threshold1.rds")
groups <- data$groups
data <- data$absolute[,fp_idx,drop=F]

scaled_counts_A <- scaled_counts_ALDEx2(data$relative, pseudocount = 0.5)
scaled_counts_D <- scaled_counts_DESeq2(data$relative, data$groups, pseudocount = 0.5)
scaled_counts_S <- scaled_counts_scran(data$relative, data$groups, pseudocount = 0.5)

calls_A <- readRDS("output/real_data_calls/no_norm/calls_oracle_ALDEx2_Owens_threshold1.rds")
calls_D <- readRDS("output/real_data_calls/no_norm/calls_oracle_DESeq2_Owens_threshold1.rds")
calls_S <- readRDS("output/real_data_calls/no_norm/calls_oracle_scran_Owens_threshold1.rds")

fp_idx <- sample(which(calls_A$rates$FP_calls & calls_D$rates$FP_calls & calls_S$rates$FP_calls), size = 1)

oracle_calls <- call_DA_NB(data$absolute[,fp_idx,drop=F], data$groups)$pval

plot_df <- data.frame(x = 1:nrow(data$relative),
                      y = data$relative[,fp_idx],
                      type = "relative abundance",
                      type2 = data$groups)
plot_df <- rbind(plot_df,
                 data.frame(x = 1:nrow(data$absolute),
                            y = data$absolute[,fp_idx],
                            type = "absolute abundance",
                            type2 = data$groups))
plot_df <- rbind(plot_df,
                 data.frame(x = 1:nrow(scaled_counts_A),
                            y = scaled_counts_A[,fp_idx],
                            type = "ALDEx2",
                            type2 = data$groups))
plot_df <- rbind(plot_df,
                 data.frame(x = 1:nrow(scaled_counts_D),
                            y = scaled_counts_D[,fp_idx],
                            type = "DESeq2",
                            type2 = data$groups))
plot_df <- rbind(plot_df,
                 data.frame(x = 1:nrow(scaled_counts_S),
                            y = scaled_counts_S[,fp_idx],
                            type = "scran",
                            type2 = data$groups))

ggplot(plot_df, aes(x = x, y = y, fill = type2)) +
  geom_point(size = 2, shape = 21) +
  facet_wrap(. ~ type, scales = "free_y") +
  theme_bw()




# Render some "misses" (discrepant genes between calls on absolute abundances
# and calls on relative abundances) while we've got the data set parsed...
# if(dataset_name == "Monaco" | dataset_name == "Klein") {
#   if(dataset_name == "Monaco") {
#     method <- "ALDEx2"
#     gene_idx <- 1228 # an arbitrary false negative (Ensa), pulled from
#     # calls_oracle_(method)_Monaco_threshold1.rds
#     scaled_counts <- scaled_counts_ALDEx2(rbind(counts_A, counts_B), pseudocount = 0.5)
#   } else {
#     method <- "scran"
#     gene_idx <- 552 # an arbitrary false positive (Cox7c), pulled from
#     # calls_oracle_(method)_Klein_threshold1.rds
#     scaled_counts <- scaled_counts_scran(rbind(counts_A, counts_B), groups, pseudocount = 0.5)
#     scaled_counts <- scaled_counts_ALDEx2(rbind(counts_A, counts_B), pseudocount = 0.5)
# 
# 
#     method <- "ALDEx2"
#     method <- "DESeq2"
#     method <- "scran"
#     gene_idx <- sample(which(readRDS(file.path("output",
#                                                "real_data_calls",
#                                                "no_norm",
#                                                paste0("calls_oracle_",method,"_",dataset_name,"_threshold1.rds")))$rates$FP_calls), size = 1)
#     scaled_counts <- scaled_counts_ALDEx2(rbind(counts_A, counts_B), pseudocount = 0.5)
#     scaled_counts <- scaled_counts_DESeq2(rbind(counts_A, counts_B), groups, pseudocount = 0.5)
#     scaled_counts <- scaled_counts_scran(rbind(counts_A, counts_B), groups, pseudocount = 0.5)
#   }
# 
#   scaled_counts_A <- scaled_counts[which(groups == levels(groups)[1]),]
#   scaled_counts_B <- scaled_counts[which(groups == levels(groups)[2]),]
# 
#   plot_df <- data.frame(x = 1:nrow(counts_A_abs),
#                         y = counts_A_abs[,gene_idx],
#                         type = "absolute abundance",
#                         type2 = levels(groups)[1])
#   plot_df <- rbind(plot_df,
#                    data.frame(x = (nrow(counts_A_abs) + 1):(nrow(counts_A_abs) + nrow(counts_B_abs)),
#                               y = counts_B_abs[,gene_idx],
#                               type = "absolute abundance",
#                               type2 = levels(groups)[2]))
#   plot_df <- rbind(plot_df,
#                    data.frame(x = 1:nrow(scaled_counts_A),
#                               # y = scaled_counts_A[,gene_idx],
#                               y = counts_A[,gene_idx],
#                               type = "scaled abundance",
#                               type2 = levels(groups)[1]))
#   plot_df <- rbind(plot_df,
#                    data.frame(x = (nrow(scaled_counts_A) + 1):(nrow(scaled_counts_A) + nrow(scaled_counts_B)),
#                               # y = scaled_counts_B[,gene_idx],
#                               y = counts_B[,gene_idx],
#                               type = "scaled abundance",
#                               type2 = levels(groups)[2]))
# 
#   if(dataset_name == "Monaco") {
#     plot_df$type2[plot_df$type2 == "CD4_naive"] <- "CD4 naive"
#   } else {
#     plot_df$type2[plot_df$type2 == "unstimulated"] <- "Unstimulated"
#   }
# 
#   temp <- data.frame(samples = 1:nrow(ref_data),
#                      group = groups,
#                      totals = rowSums(ref_data),
#                      type = "absolute")
#   temp <- rbind(temp,
#                 data.frame(samples = 1:nrow(data),
#                            group = groups,
#                            totals = rowSums(data),
#                            type = "relative"))
#   temp <- rbind(temp,
#                 data.frame(samples = 1:(nrow(scaled_counts_A) + nrow(scaled_counts_B)),
#                            group = groups,
#                            totals = c(rowSums(scaled_counts_A), rowSums(scaled_counts_B)),
#                            type = "rescaled"))
#   ggplot(temp, aes(x = samples, y = totals, fill = group)) +
#     geom_bar(stat = "identity") +
#     facet_wrap(. ~ type, scales = "free_y")
# 
#   p <- ggplot(plot_df, aes(x = x, y = y, fill = type2)) +
#     geom_point(size = 2, shape = 21) +
#     facet_wrap(. ~ type) +
#     theme_bw() +
#     labs(x = "sample index",
#          y = "abundance",
#          fill = "Condition or type") +
#     scale_fill_manual(values = c("#0771DE", "orange"))
# 
#   show(p)
# 
#   if(dataset_name == "Monaco") {
#     ggsave("output/images/FN_Monaco.png",
#            p,
#            dpi = 100,
#            units = "in",
#            height = 2.5,
#            width = 6)
#   } else if(dataset_name == "Klein") {
#     show(p)
#     ggsave("output/images/FP_Klein.png",
#            p,
#            dpi = 100,
#            units = "in",
#            height = 2.5,
#            width = 6)
#   }
# }
