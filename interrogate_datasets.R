source("path_fix.R")

library(RSQLite)
library(ggplot2)

library(optparse)

option_list = list(
  make_option(c("--p"), type = "numeric", default = 100,
              help = "number of features", metavar = "numeric"),
  make_option(c("--corrp"), type = "numeric", default = 0,
              help = "correlation flag", metavar = "numeric")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

p <- opt$p
corrp <- opt$corrp

# Testing
# p <- 100
# corrp <- 0

# Pull UUIDs matching P, CORRP characterstics
conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))
res <- dbGetQuery(conn, paste0("SELECT UUID FROM datasets WHERE P=",p," AND CORRP=",corrp))
fc_distro <- c()
de_distro <- c()
for(i in 1:length(res$UUID)) {
  cat(paste0("Parsing dataset ", i, " / ", length(res$UUID), "\n"))
  uuid <- res$UUID[i]
  data <- readRDS(file.path("output", "datasets", paste0(uuid, ".rds")))
  # Calculate FC
  M <- data$simulation$abundances
  counts_A <- M[1:(nrow(M)/2),]
  counts_B <- M[(nrow(M)/2+1):nrow(M),]
  m1 <- mean(rowSums(counts_A))
  m2 <- mean(rowSums(counts_B))
  fc <- max(c(m1, m2)) / min(c(m1, m2))
  fc_distro <- c(fc_distro, fc)
  # Estimate % DE features
  m1 <- colMeans(counts_A + 0.1)
  m2 <- colMeans(counts_B + 0.1)
  fc <- m1 / m2
  prop_de <- sum(fc <= (2/3) | fc >= (3/2)) / length(fc)
  de_distro <- c(de_distro, prop_de)
}

plot_df <- data.frame(x = fc_distro)
pl <- ggplot(plot_df, aes(x = x)) +
  geom_histogram(color = "#ffffff") +
  labs(x = "fold change across conditions")
ggsave(paste0("data_hist_fc_p",p,"_corrp",corrp,".png"), pl, units = "in", dpi = 100, height = 4, width = 5)

plot_df <- data.frame(x = de_distro)
pl <- ggplot(plot_df, aes(x = x)) +
  geom_histogram(color = "#ffffff") +
  labs(x = "proportion differentially abundant features")
ggsave(paste0("data_hist_de_p",p,"_corrp",corrp,".png"), pl, units = "in", dpi = 100, height = 4, width = 5)

