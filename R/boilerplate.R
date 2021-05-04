#' Generate a color palette suitable for stacked bar plots given a feature number (S)
#'
#' @param S number of features (colors to generate)
#' @return list of hex colors
#' @import RColorBrewer
#' @export
generate_highcontrast_palette <- function(S) {
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  sample(getPalette(S))
}

#' Generate a color palette suitable for stacked bar plots given a feature number (S)
#'
#' @param data matrix of abundances or observed counts (samples x features)
#' @param palette a color palette
#' @param save_name optional name to save plot under
#' @import RColorBrewer
#' @export
plot_stacked_bars <- function(data, palette = NULL, save_name = NULL) {
  data <- as.data.frame(t(data))
  # data is now taxa x samples
  n_samples <- ncol(data)
  data <- cbind(1:nrow(data), data)
  colnames(data) <- c("feature", 1:n_samples)
  data_long <- pivot_longer(data, !feature, names_to = "sample", values_to = "abundance")
  data_long$feature <- factor(data_long$feature)
  data_long$sample <- factor(data_long$sample, levels = 1:n_samples)
  
  if(is.null(palette)) {
    palette <- generate_highcontrast_palette(length(unique(data_long$feature)))
  }
  
  p <- ggplot(data_long, aes(fill = feature, y = abundance, x = sample)) + 
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = palette) +
    theme_bw() +
    theme(legend.position = "none")
  show(p)
  if(!is.null(save_name)) {
    ggsave(file.path("output", "images", save_name),
           p,
           units = "in",
           dpi = 100,
           height = 5,
           width = 6)
  }
}

#' Generate a bipartite graph for differential abundance/expression
#'
#' @param expr1 expression vector in baseline condition
#' @param expr2 expression vector in differential condition
#' @param alpha optional transparency parameter
#' @return NULL
#' @import ggplot2
#' @export
plot_bipartite_graph <- function(expr1, expr2, alpha = 1) {
  # Simulate data
  # n <- 100
  # data <- data.frame(x = 1, xend = 2, y = runif(n, min = 0, max = 100))
  # data$yend <- rnorm(n, data$y, 20)
  
  plot_data <- data.frame(x = 1, xend = 2, y = expr1, yend = expr2)
  qq <- quantile(plot_data$y, probs = seq(from = 0, to = 1, length.out = 5))
  qq[[1]] <- -Inf
  qq[[length(qq)]] <- Inf
  
  plot_data$rank <- cut(plot_data$y, qq) # factor
  
  p <- ggplot(plot_data, aes(color = rank)) +
    geom_point(aes(x = x, y = y), size = 3, alpha = alpha) +
    geom_point(aes(x = xend, y = yend), size = 3, alpha = alpha) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend), alpha = alpha) +
    scale_x_discrete(name = "condition", 
                     limits = c("baseline", "differential")) +
    ylab("abundance") +
    theme_bw() +
    theme(legend.position = "none")
  show(p)
}
