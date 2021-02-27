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
