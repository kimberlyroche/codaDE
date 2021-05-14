source("path_fix.R")

library(stringr)

file.dirs <- c("p100_corrp0", "p100_corrp0.5", "p100_corrp1",
               "p1000_corrp0", "p1000_corrp0.5", "p1000_corrp1",
               "p5000_corrp0", "p5000_corrp0.5", "p5000_corrp1")

for(file.dir in file.dirs) {
  cat("Parsing files in:", file.dir, "\n")
  # Generalize to take a `p` argument
  # files <- list.files(pattern = "simresults_p100_simulated_*", path = "output")
  files <- list.files(path = file.path("output", file.dir), pattern = paste0("simresults_.*?_\\w{8}-.*"))
  results_all <- NULL
  for(file in files) {
    results <- readRDS(file.path("output", file.dir, file))
    if(is.null(results_all)) {
      results_all <- results
    } else {
      results_all <- rbind(results_all, results)
    }
  }
  saveRDS(results_all, file = file.path("output", file.dir, paste0("simresults_", file.dir, "_all.rds")))
}
