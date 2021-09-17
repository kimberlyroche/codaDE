data <- readRDS(file.path("output", "calls_summary.rds"))

# data <- data %>%
#   filter((dataset == "Monaco" & threshold == 2) | (dataset != "Monaco" & threshold == 1))

for(i in 1:5) {
  cat(paste0("  TPR: ", round(mean(data %>% filter(threshold == i) %>% pull(TPR)), 2), "\n"))
  cat(paste0("1-FPR: ", round(mean(1 - data %>% filter(threshold == i) %>% pull(FPR)), 2), "\n\n"))
}

data$dataset <- factor(data$dataset,
                       levels = c("VieiraSilva",
                                  "Barlow",
                                  "Song",
                                  "Monaco",
                                  "Hagai",
                                  "Owens",
                                  "Klein",
                                  "Yu"))

data %>%
  mutate(ds_id = as.numeric(dataset)) %>%
  mutate(TPR = round(TPR, 2)) %>%
  mutate(FPR = 1 - round(FPR, 2)) %>%
  arrange(method, ds_id)

ggplot(data, aes(x = TPR, color = method)) +
  geom_histogram() +
  facet_wrap(. ~ method)

ggplot(data, aes(x = 1 - FPR, color = method)) +
  geom_histogram() +
  facet_wrap(. ~ method)

data %>%
  mutate(TPR = round(TPR, 2)) %>%
  mutate(FPR = 1 - round(FPR, 2)) %>%
  group_by(dataset) %>%
  summarize(mTPR = mean(TPR))

data %>%
  mutate(TPR = round(TPR, 2)) %>%
  mutate(FPR = 1 - round(FPR, 2)) %>%
  group_by(dataset) %>%
  summarize(mFPR = mean(FPR))
