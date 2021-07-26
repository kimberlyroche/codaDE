source("path_fix.R")

library(codaDE)
library(RSQLite)
library(dplyr)
library(optparse)

option_list = list(
  make_option(c("--p"), type = "numeric", default = 100,
              help = "number of genes", metavar = "numeric")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

p <- opt$p

corrp_sweep <- c(0, 1, 2, 3, 4)
log_mean_sweep <- c(2, 3, 4, 5, 6)
perturbation_sweep <- c(0.2, 0.4, 0.6, 0.8, 1) # percent log mean
# perturbation_sweep <- c(0.5, 1, 2, 3, 4)
# perturbation_sweep <- c(0.1, 0.2, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)
rep_noise_sweep <- c(0.2, 0.6, 1)
# rep_noise_sweep <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
diff_sweep <- c(0.1, 0.3, 0.5, 0.7)

n_combos <- length(corrp_sweep)*length(log_mean_sweep)*length(perturbation_sweep)*length(rep_noise_sweep)*length(diff_sweep)

# Create DB (if doesn't exist?)
conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

results <- dbGetQuery(conn, "SELECT * FROM datasets")
settings <- data.frame(ID = c(),
                       P = c(),
                       CORRP = c(),
                       LOG_MEAN = c(),
                       PERTURBATION = c(),
                       REP_NOISE = c(),
                       PERCENT_DIFF = c())

counter <- 1
for(corrp in corrp_sweep) {
  for(log_mean in log_mean_sweep) {
    for(perturbation in perturbation_sweep) {
      for(replicate_noise in rep_noise_sweep) {
        for(percent_diff in diff_sweep) {
          if(counter %% 100 == 0) {
            cat(paste0("Iteration ", counter, " / ", n_combos, "\n"))
          }
          counter <- counter + 1
          uuids <- results %>%
            filter(P == p & CORRP == corrp &
                     LOG_MEAN == log_mean &
                     PERTURBATION == perturbation &
                     REP_NOISE == replicate_noise &
                     PERCENT_DIFF_SIM == percent_diff) %>%
            pull(UUID)
          if(length(uuids) > 0) {
            next
          }
          # Add this setting
          settings <- rbind(settings,
                            data.frame(ID = counter,
                                       P = p,
                                       CORRP = corrp,
                                       LOG_MEAN = log_mean,
                                       PERTURBATION = perturbation,
                                       REP_NOISE = replicate_noise,
                                       PERCENT_DIFF = percent_diff))
        }
      }
    }
  }
}

write.table(settings, file = paste0("input_gen_", p, ".txt"), sep = "\t")

dbDisconnect(conn)
