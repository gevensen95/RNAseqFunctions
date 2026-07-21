#' Ol' Reliable ggplot2 theme
#'
#' A consistent ggplot2 theme used across RNAseqBasic plots: clean panel
#' borders, subtle gridlines, and bold black-and-white facet strips.
#'
#' @return A ggplot2 theme object to be added to a ggplot.
#' @export
Ol_Reliable <- function() {
  ggplot2::theme(
    panel.border = ggplot2::element_rect(colour = "black", fill = NA),
    panel.background = ggplot2::element_rect(fill = "white"),
    panel.grid.major.x = ggplot2::element_line(colour = "grey", linewidth = 0.25),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(colour = "grey", linewidth = 0.25),
    panel.grid.minor.y = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(colour = "black"),
    axis.title = ggplot2::element_text(colour = "black", face = "italic"),
    axis.ticks = ggplot2::element_line(colour = "black"),
    legend.title = ggplot2::element_text(hjust = 0.5),
    strip.background = ggplot2::element_rect(fill = "black"),
    strip.text = ggplot2::element_text(colour = "white", face = "bold"),
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = ggplot2::element_text(hjust = 0.5),
    legend.key = ggplot2::element_rect(fill = "white")
  )
}
