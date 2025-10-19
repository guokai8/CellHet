#' CellHet Default Color Palette
#'
#' Returns the default color palette used across all CellHet visualizations
#'
#' @param n Number of colors to return. If NULL, returns all colors.
#' @param alpha Transparency level (0-1)
#'
#' @return A character vector of hex color codes
#' @export
#'
#' @examples
#' # Get all colors
#' cellhet_colors()
#'
#' # Get first 10 colors
#' cellhet_colors(10)
#'
#' # Get colors with transparency
#' cellhet_colors(5, alpha = 0.7)
cellhet_colors <- function(n = NULL, alpha = 1) {

  # Default CellHet color palette
  colors <- c(
    "#2670B8", "#B3D9CE", "#16A16E",
    "#876AAA", "#A5CDA7", "#9E1B8E", "#B85F81",
    "#405BA0", "#C0F5D7", "#223C4D", "#4758A0", "#3F578B", "#C2DF7E", "#71B9A3",
    "#BDFBDA", "#D7E7D1", "#F8C4C9", "#E79292", "#6C9D8E", "#58B390", "#B6DDCB", "#E0F0E5",
    "#16416E", "#5499B2", "#83BBB7", "#B5D6B9", "#DAE6D3", "#6FB1B9", "#A8D7D5", "#B7CDD1",
    "#82BCD8", "#AFE8F3", "#B5F8F1", "#494574", "#427794", "#4DA8A0", "#94CFB2", "#6D8C87",
    "#87B6AA", "#C9E1BE", "#4D85A0", "#82BED8", "#C2DFE7", "#91BE93", "#BFD8D4", "#D35400",
    "#E67E22", "#F1C40F", "#FBC02D", "#EC407A", "#F06292", "#00ACC1", "#26C6DA", "#C0392B"
  )

  # Apply transparency if needed
  if (alpha < 1) {
    colors <- adjustcolor(colors, alpha.f = alpha)
  }

  # Return requested number of colors
  if (is.null(n)) {
    return(colors)
  } else {
    # Cycle through colors if more are needed
    return(rep_len(colors, n))
  }
}

#' Scale Color for CellHet
#'
#' ggplot2 color scale using CellHet default palette
#'
#' @param ... Additional arguments passed to scale_color_manual
#'
#' @return A ggplot2 scale object
#' @export
scale_color_cellhet <- function(...) {
  ggplot2::scale_color_manual(values = cellhet_colors(), ...)
}

#' Scale Fill for CellHet
#'
#' ggplot2 fill scale using CellHet default palette
#'
#' @param ... Additional arguments passed to scale_fill_manual
#'
#' @return A ggplot2 scale object
#' @export
scale_fill_cellhet <- function(...) {
  ggplot2::scale_fill_manual(values = cellhet_colors(), ...)
}
