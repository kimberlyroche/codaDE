# Parse simulated data and analysis and write summary results to the file "results.tsv"

# use only this number of random outputs
# limit <- 1000

simulations <- list.files(path = "simulated_data", pattern = "*.rds", full.names = TRUE, recursive = FALSE)
# simulations <- simulations[sample(1:length(simulations))[1:limit]]

out_file <- file("results.txt", 'w')
cat(paste("p",
            "prop_de",
            "sf_corr",
            "expr_de",
            "dir_de",
            "sparsity",
            "entropy",
            "tp",
            "tn",
            "fp",
            "fn", sep = "\t"), "\n", file = out_file)

for(i in 1:length(simulations)) {
  if(i %% 100 == 0) {
    cat("Parsing simulation #",i,"\n")
  }
  simdat <- readRDS(simulations[i])$result
  if(!is.null(simdat)) {
    cat(paste(simdat$p,
                simdat$proportion_da,
                simdat$size_factor_correlation,
                simdat$median_expr_da_quantile,
                simdat$net_dir_da,
                simdat$sparsity,
                simdat$sim_entropy,
                simdat$tp,
                simdat$tn,
                simdat$fp,
                simdat$fn, sep = "\t"), "\n", file = out_file, append = TRUE)
  }
}

close(out_file)

