library(RSQLite)
# library(GauPro)
library(mlegp)
library(tidyverse)
library(gridExtra)

# ------------------------------------------------------------------------------
#   Testing DB setup and query
# ------------------------------------------------------------------------------

# Load and clean up cars data set
data("mtcars")
mtcars$car_names <- rownames(mtcars)
rownames(mtcars) <- c()
head(mtcars)

# Create DB (if doesn't exist?)
conn <- dbConnect(RSQLite::SQLite(), "temp.db")

# List tables
dbListTables(conn)

# Write full data set to new table
dbWriteTable(conn, "cars_data", mtcars)
dbListTables(conn)

dbGetQuery(conn, "SELECT mpg FROM cars_data")
dbExecute(conn, "INSERT into cars_data VALUES(21.0,6,200,100,3.5,2.5,15.00,1,0,4,1,'Fake Car')")
dbGetQuery(conn, "SELECT * FROM cars_data WHERE car_names = 'Fake Car'")

dbDisconnect(conn)

# ------------------------------------------------------------------------------
#   Testing GP regression on noise data
#
#   We want multidimensional input to multidimensional output. The mlegp package
#   seems to do that.
# ------------------------------------------------------------------------------

plot_2x2 <- function(data_true, data_pred = NULL) {
  p1 <- ggplot()
  if(!is.null(data_pred)) {
    p1 <- p1 + geom_ribbon(data = data_pred, mapping = aes(x = x,
                                                           ymin = z1 - 2*z1_se,
                                                           ymax = z1 + 2*z1_se),
                           color = "gray",
                           alpha = 0.5) +
      geom_line(data = data_pred, mapping = aes(x = x, y = z1))
  }
  p1 <- p1 + geom_point(data = data_true, mapping = aes(x = x, y = z1))
  p2 <- ggplot()
  if(!is.null(data_pred)) {
    p2 <- p2 + geom_ribbon(data = data_pred, mapping = aes(x = x,
                                                           ymin = z2 - 2*z2_se,
                                                           ymax = z2 + 2*z2_se),
                           color = "gray",
                           alpha = 0.5) +
      geom_line(data = data_pred, mapping = aes(x = x, y = z2))
  }
  p2 <- p2 + geom_point(data = data_true, mapping = aes(x = x, y = z2))
  p3 <- ggplot()
  if(!is.null(data_pred)) {
    p3 <- p3 + geom_ribbon(data = data_pred, mapping = aes(x = y,
                                                           ymin = z1 - 2*z1_se,
                                                           ymax = z1 + 2*z1_se),
                           color = "gray",
                           alpha = 0.5) +
      geom_line(data = data_pred, mapping = aes(x = y, y = z1))
  }
  p3 <- p3 + geom_point(data = data_true, mapping = aes(x = y, y = z1))
  p4 <- ggplot()
  if(!is.null(data_pred)) {
    p4 <- p4 + geom_ribbon(data = data_pred, mapping = aes(x = y,
                                                           ymin = z2 - 2*z2_se,
                                                           ymax = z2 + 2*z2_se),
                           color = "gray",
                           alpha = 0.5) +
      geom_line(data = data_pred, mapping = aes(x = y, y = z2))
  }
  p4 <- p4 + geom_point(data = data_true, mapping = aes(x = y, y = z2))
  grid.arrange(grobs = list(p1, p2, p3, p4), ncol = 2)
}

# ------------------------------------------------------------------------------
#   Simulate data
# ------------------------------------------------------------------------------

# 2D -> 2D function
f <- function(x, y) {
  list(z1 = x**2 + y, z2 = 2*y - sqrt(x))
}

input_true <- data.frame(x = seq(from = 0, to = 2, length.out = 20),
                         y = seq(from = -1, to = 1, length.out = 20))

output_true <- as.data.frame(f(input_true$x, input_true$y))
# Add noise
output_true$z1 <- output_true$z1 + rnorm(nrow(output_true), 0, 0.5)
output_true$z2 <- output_true$z2 + rnorm(nrow(output_true), 0, 0.5)

data_true <- bind_cols(input_true, output_true)
plot_2x2(data_true)

# ------------------------------------------------------------------------------
#   Fit and predict 2D -> 2D GP
# ------------------------------------------------------------------------------

gp <- mlegp(input_true, output_true)

# Interpolate
p_len <- 100
input_pred <- data.frame(x = seq(from = 0, to = 2, length.out = p_len),
                         y = seq(from = -1, to = 1, length.out = p_len))
output_pred_z1 <- predict(gp[[1]], newData = input_pred, se.fit = TRUE)
output_pred_z2 <- predict(gp[[2]], newData = input_pred, se.fit = TRUE)
output_pred <- data.frame(z1 = output_pred_z1$fit,
                          z2 = output_pred_z2$fit,
                          z1_se = output_pred_z1$se.fit,
                          z2_se = output_pred_z2$se.fit)

data_pred <- bind_cols(input_pred, output_pred)
plot_2x2(data_true, data_pred)
