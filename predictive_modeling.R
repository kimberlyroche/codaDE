source("path_fix.R")

library(codaDE)
library(RSQLite)

# Create DB (if doesn't exist?)
conn <- dbConnect(RSQLite::SQLite(), file.path("output", "simulations.db"))

# Pull data characterstics and results (join)

# Define test/train set

# Fit GP on train

# Evaluate on test (visualize)

dbDisconnect(conn)
