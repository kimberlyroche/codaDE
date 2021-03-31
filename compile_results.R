if(R.version$major == 4) {
  cat("Updating lib paths...\n")
  .libPaths(c("/gpfs/fs1/data/mukherjeelab/roche/Rlibs", .libPaths()[2]))
}

library(stringr)

file.dir <- "p15000"

# Generalize to take a `p` argument
# files <- list.files(pattern = "simresults_p5000_simulated_*", path = "output")
files <- list.files(path = file.path("output", file.dir), pattern = paste0("simresults_", file.dir, "_simulated_bulk_\\w{8}-.*"))
results_all <- NULL
for(file in files) {
  results <- readRDS(file.path("output", file.dir, file))
  if(is.null(results_all)) {
    results_all <- results
  } else {
    results_all <- rbind(results_all, results)
  }
}
saveRDS(results_all, file = file.path("output", file.dir, paste0("simresults_", file.dir, "_simulated_all.rds")))
