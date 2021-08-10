source("path_fix.R")

library(codaDE)
library(RSQLite)
library(rfinterval)

# ------------------------------------------------------------------------------
#   Test and visualize prediction intervals for RF
# ------------------------------------------------------------------------------

n <- 1000
data <- data.frame(rnorm(n))
for(i in 1:19) {
  data <- cbind(data, rnorm(n))
}
colnames(data) <- paste0("X", 1:ncol(data))
Y <- rowSums(data)
data$Y <- Y + rnorm(n, 0, 0.01)
n_train <- round(n/2)
train_idx <- sample(1:n, size = n_train)
output <- rfinterval(Y ~ ., train_data = data[train_idx,], test_data = data[-train_idx,],
                     method = c("oob", "split-conformal", "quantreg"),
                     symmetry = TRUE, alpha = 0.1)
cat(paste0("Percent of observations covered by 90% interval: ",
           round(sum(output$oob_interval$lo < output$test_data$Y & output$oob_interval$up > output$test_data$Y) / n_train, 3)*100,
           "\n"))

# It really misses the outliers. Predictions are pretty flat overall.
plot_df <- data.frame(lower = output$oob_interval$lo,
                      true_val = output$test_data$Y,
                      upper = output$oob_interval$up)
plot_df <- plot_df %>%
  arrange(true_val)
ggplot(plot_df, aes(x = 1:nrow(plot_df), y = true_val, ymin = lower, ymax = upper)) +
  geom_ribbon(fill = "grey70") +
  geom_line()

# ------------------------------------------------------------------------------
#   Visualize partially accurate abundance simulations that capture at least
#   some of the true fold change between conditions
# ------------------------------------------------------------------------------

conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

# Cases where there was a net DECREASE in abundance
# (Note: relative abundances should always be around 1 FC)
res <- dbGetQuery(conn, "SELECT P, FC_ABSOLUTE, FC_PARTIAL, FC_RELATIVE FROM datasets WHERE FC_ABSOLUTE < 1")
plot_df <- data.frame(rel = res$FC_RELATIVE,
                      partial = res$FC_PARTIAL,
                      abs = res$FC_ABSOLUTE)
plot_df <- plot_df %>%
  filter(abs > 0.5) %>%
  filter(abs < partial & partial < rel) %>%
  arrange(abs)
ggplot(plot_df, aes(x = nrow(plot_df):1, ymin = abs, y = partial, ymax = rel)) +
  geom_ribbon(fill = "grey70") +
  geom_line()

cat(paste0(nrow(res), " decrease cases\n"))

# Cases where there was a net INCREASE in abundance
res <- dbGetQuery(conn, "SELECT P, FC_ABSOLUTE, FC_PARTIAL, FC_RELATIVE FROM datasets WHERE FC_ABSOLUTE > 1")
plot_df <- data.frame(rel = res$FC_RELATIVE,
                      partial = res$FC_PARTIAL,
                      abs = res$FC_ABSOLUTE)
plot_df <- plot_df %>%
  filter(abs < 2) %>%
  filter(rel < partial & partial < abs) %>%
  arrange(abs)
ggplot(plot_df, aes(x = 1:nrow(plot_df), ymin = abs, y = partial, ymax = rel)) +
  geom_ribbon(fill = "grey70") +
  geom_line()

cat(paste0(nrow(res), " increase cases\n"))

dbDisconnect(conn)
