library(dplyr)

md <- read.table(file.path("simulated_data", "metadata.tsv"), header = TRUE, stringsAsFactors = FALSE)
filter_abundance <- 10
p <- 10000

# (1) Show all recorded runs associated with all conditions for a given gene number.
temp <- md[md$filter_threshold == filter_abundance & md$p == p,] %>%
  group_by(p, proportion_da, size_factor_correlation) %>%
  tally()
print(as_tibble(temp), n = 100)

# (2) Get conditions from a set of runs with more or less than 20 replicates.
#     This should be empty!
temp <- md[md$filter_threshold == filter_abundance,] %>%
  group_by(p, proportion_da, size_factor_correlation) %>%
  tally() %>%
  filter(n != 20)
temp

# If it's not empty, dump out argument strings associated with missing conditions; we can feed these
# to an instance of `run.R` to fill in missing instances.
use_existing_str <- "TRUE"
if(nrow(temp) > 0) {
  temp <- as.data.frame(temp)
  temp <- temp[order(temp$n),]
  for(i in 1:nrow(temp)) {
    cat(paste0("--p=",temp[i,]$p," --n=250 --k=1 --prop_da=",temp[i,]$proportion_da," --sf_corr=",temp[i,]$size_factor_correlation," --filter_abundance=",filter_abundance," --existing=",use_existing_str,"  (",20 - temp[i,]$n,")\n"))
  }
}

