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
#' @param data a data.frame with (named) columns for feature, sample, and abundance
#' @param palette a color palette
#' @param save_name optional name to save plot under
#' @import RColorBrewer
#' @export
plot_stacked_bars <- function(data, palette = NULL, save_name = NULL) {
  if(is.null(palette)) {
    palette <- generate_highcontrast_palette(length(unique(data$feature)))
  }
  p <- ggplot(data, aes(fill = feature, y = abundance, x = sample)) + 
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = palette) +
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

