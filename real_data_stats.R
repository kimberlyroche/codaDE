library(tidyverse)

datasets <- c("VieiraSilva", "Barlow", "Song", "Monaco", "Hagai", "Owens", "Klein", "Yu")
thresholds <- c(1, 1, 1, 1, 1, 1, 1, 1)

bar_components <- NULL
compiled_stats <- NULL
accuracy <- c()
fpr <- c()
tpr <- c()
for(i in 1:length(datasets)) {
  for(method in c("ALDEx2", "DESeq2", "scran")) {
    if(datasets[i] == "Monaco" & method == "scran") next;
    calls <- readRDS(file.path("output",
                               "real_data_calls",
                               "no_norm",
                               paste0("calls_oracle_",method,"_",datasets[i],"_threshold",thresholds[i],".rds")))
    TP <- sum(calls$rates$TP_calls)
    TN <- sum(calls$rates$TN_calls)
    FP <- sum(calls$rates$FP_calls)
    FN <- sum(calls$rates$FN_calls)
    
    bar_components <- rbind(bar_components,
                            data.frame(dataset = datasets[i], method = method, count = TP, type = "TP"))
    bar_components <- rbind(bar_components,
                            data.frame(dataset = datasets[i], method = method, count = FN, type = "FN"))
    bar_components <- rbind(bar_components,
                            data.frame(dataset = datasets[i], method = method, count = FP, type = "FP"))
    bar_components <- rbind(bar_components,
                            data.frame(dataset = datasets[i], method = method, count = TN, type = "TN"))
    
    accuracy <- c(accuracy,
                  (TP+TN)/(TP+TN+FP+FN))
    fpr <- c(fpr, FP/(FP+TN))
    tpr <- c(tpr, TP/(TP+FN))
    compiled_stats <- rbind(compiled_stats,
                            data.frame(dataset = datasets[i],
                                       method = method,
                                       sensitivity = TP/(TP+FN),
                                       specificity = 1 - FP/(FP+TN)))
  }
}

bar_components$dataset <- factor(bar_components$dataset, levels = c("VieiraSilva",
                                                                    "Barlow",
                                                                    "Song",
                                                                    "Monaco",
                                                                    "Hagai",
                                                                    "Owens",
                                                                    "Klein",
                                                                    "Yu"))
levels(bar_components$dataset) <- c("Vieira-Silva et al.",
                                    "Barlow et al.",
                                    "Song et al.",
                                    "Monaco et al.",
                                    "Hagai et al.",
                                    "Owens et al.",
                                    "Klein et al.",
                                    "Yu et al.")
bar_components$type <- factor(bar_components$type, levels = c("TP", "TN", "FN", "FP"))
levels(bar_components$type) <- c("True positive", "True negative", "False negative", "False positive")

# Plot barplot
p <- ggplot(bar_components, aes(x = method, y = count, fill = type, label = count)) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  facet_wrap(. ~ dataset, scales = "free", ncol = 4) +
  theme_bw() +
  scale_fill_brewer(palette = "GnBu") +
  labs(fill = "",
       x = "Differential abundance calling method")

ggsave("output/temp2.png",
       p,
       dpi = 100,
       units = "in",
       height = 6,
       width = 12)

min(compiled_stats$sensitivity)
median(compiled_stats$sensitivity)
max(compiled_stats$sensitivity)

min(compiled_stats %>% filter(method != "ALDEx2") %>% pull(sensitivity))
median(compiled_stats %>% filter(method != "ALDEx2") %>% pull(sensitivity))
max(compiled_stats %>% filter(method != "ALDEx2") %>% pull(sensitivity))

compiled_stats %>% group_by(dataset) %>% summarize(ms = mean(specificity)) %>% arrange(desc(ms))


