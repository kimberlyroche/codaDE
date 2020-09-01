simdata_files <- list.files(path = "simulated_data", pattern = "*.rds")

args = commandArgs(trailingOnly=TRUE)
lower <- as.numeric(args[1])
upper <- as.numeric(args[2])

for(i in lower:upper) {
  filename <- simdata_files[i]
  data <- tryCatch({
    readRDS(file.path("simulated_data", filename))
  },
  warning = function(w) { cat("Warning -- could not open file:",filename,"\n") },
  error = function(e) { cat("Error -- could not open file:",filename,"\n") }
  )
  if(!is.null(data)) {
    cat(filename,"\t",data$properties$p,"\t",data$properties$proportion_da,"\t",data$properties$size_factor_correlation,"\n")
    res1 <- data$results[[1]]
    res2 <- data$results[[2]]
    saveRDS(res1, file = file.path("simulated_analyses", "analysis1", filename))
    if(res2$filter_abundance == 10) {
      saveRDS(res2, file = file.path("simulated_analyses", "analysis2", filename))
    } else if(res2$filter_abundance == 3) {
      saveRDS(res2, file = file.path("simulated_analyses", "analysis3", filename))
    }
  }
}

