source("path_fix.R")

library(stringr)

results <- NULL

file_list <- list.files(path = "output",
                        pattern = "calls_self_.*\\.rds",
                        full.names = TRUE)
for(file in file_list) {
  data <- readRDS(file)
  str_pieces <- str_match(file, "calls_self_(.*?)_(.*?)_threshold(\\d+)\\.rds")
  results <- rbind(results,
                   data.frame(dataset = str_pieces[1,3],
                              method = str_pieces[1,2],
                              threshold = as.numeric(str_pieces[1,4]),
                              TPR = data$rates$TPR,
                              FPR = data$rates$FPR))
}

saveRDS(results, file.path("output", "calls_summary.rds"))
