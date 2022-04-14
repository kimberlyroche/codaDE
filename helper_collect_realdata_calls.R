source("path_fix.R")

library(stringr)

results <- NULL
file_list <- list.files(path = file.path("output", "real_data_calls", "no_norm"),
                        pattern = "calls_oracle_.*\\.rds",
                        full.names = TRUE)
for(file in file_list) {
  data <- readRDS(file)
  str_pieces <- str_match(file, "calls_oracle_(.*?)_(.*?)_threshold1\\.rds")
  results <- rbind(results,
                   data.frame(dataset = str_pieces[1,3],
                              method = str_pieces[1,2],
                              TPR = data$rates$TPR,
                              FPR = data$rates$FPR))
}

results$TPR <- round(results$TPR, 3)
results$FPR <- round(results$FPR, 3)

results %>% filter(dataset == "VieiraSilva")
results %>% filter(dataset == "Muraro")
results %>% filter(dataset == "Hagai")
results %>% filter(dataset == "Hashimshony")
results %>% filter(dataset == "Gruen")
results %>% filter(dataset == "Kimmerling")

results %>% filter(dataset == "Song")
results %>% filter(dataset == "Barlow")
results %>% filter(dataset == "Monaco")
results %>% filter(dataset == "Yu")
results %>% filter(dataset == "Klein")
results %>% filter(dataset == "Owens")

saveRDS(results, file.path("output", "calls_summary.rds"))
