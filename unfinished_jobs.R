library(dplyr)

# V1: simulated_data directory
md <- read.table(file.path("simulated_data", "metadata.tsv"), header = TRUE, stringsAsFactors = FALSE)

# (1) Show all recorded runs associated with all conditions for a GIVEN GENE NUMBER
temp <- md %>%
  group_by(p, proportion_da, size_factor_correlation) %>%
  tally()
print(as_tibble(temp), n = 100)

# (2) Get conditions from a set of runs with more or less than 20 replicates.
#     This should be empty!
temp <- md %>%
  group_by(p, proportion_da, size_factor_correlation) %>%
  tally() %>%
  filter(n != 20)
temp

# If it's not empty, dump out argument strings associated with missing conditions; we can feed these
# to an instance of `run.R` to fill in missing instances.
use_existing_str <- "TRUE"
new_filter_abundance <- 3
save_slot <- 2
p <- 10000
if(nrow(temp) > 0) {
  temp <- as.data.frame(temp)
  temp <- temp[order(temp$n),]
  for(i in 1:nrow(temp)) {
    if(temp[i,]$p == p) {
      arguments_str <- paste0("--p=",temp[i,]$p," --n=250 --k=1 --prop_da=",temp[i,]$proportion_da," --sf_corr=",temp[i,]$size_factor_correlation," --filter_abundance=",new_filter_abundance," --existing=",use_existing_str," --save_slot=",save_slot)
      # cat(arguments_str,"\n")
      fileConn <- file("job.slurm")
      writeLines(c('#!/bin/bash',
        paste0('#SBATCH -J DA',temp[i,]$p),
        '#SBATCH --mem=16GB',
        '#SBATCH --get-user-env',
        '#SBATCH --time=00:30:00',
        '#\n',
        'module add R/3.6.1-gcb03',
        'module add gcc/7.1.0-fasrc01\n',
        'cd /data/mukherjeelab/roche/codaDE',
        paste0('srun Rscript run.R ',arguments_str)), fileConn)
      close(fileConn)
      submit_str <- paste0("sbatch --array=1-",temp[i,]$n," job.slurm")
      # cat(submit_str,"\n")
      system(submit_str)
    }
  }
}

# V2: simulated_analysis directory
analysis_label <- "analysis3"
filenames <- list.files(path = file.path("simulated_analyses", analysis_label), pattern = "*.rds")
md <- read.table(file.path("simulated_data", "metadata.tsv"), header = TRUE, stringsAsFactors = FALSE)
md <- md[md$filename %in% filenames,]
temp <- md %>%
  group_by(p, proportion_da, size_factor_correlation) %>%
  tally() %>%
  filter(n != 20)
print(as_tibble(temp), n = 100)
