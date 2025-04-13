#' Create a Gene Regulation Dotplot
#'
#' Creates a dotplot visualizing up and down regulated genes across different comparisons
#' with customizable shapes and count label options.
#'
#' @param data A data frame containing gene regulation data
#' @param comparison_col Column name for the comparison categories
#' @param cell_type_col Column name for the cell types
#' @param regulation_col Column name for the regulation direction (should contain "Up" and "Down" values)
#' @param count_col Column name for the gene counts
#' @param count_display Where to display count numbers: "none", "inside", "below", or "both" (default = "none")
#' @param use_triangles Logical. If TRUE, uses triangles (shapes 24, 25), otherwise uses circles (shape 21) (default = TRUE)
#' @param point_size_range Numeric vector of length 2 for the range of point sizes (default = c(1, 10))
#' @param up_color Color for up-regulated genes (default = "deepskyblue")
#' @param down_color Color for down-regulated genes (default = "darkorange")
#' @param offset Size of horizontal offset between up/down points (default = 0.2)
#' @param title Plot title (default = "Up and Down Regulated Genes Across Comparisons")
#' @param x_angle Angle for x-axis text (default = 45)
#'
#' @return A ggplot object
#'
#' @examples
#' # Create a sample dataset
#' sample_data <- data.frame(
#'   comparison = rep(c("Day14_WT_vs_Day7_WT", "Day7_TF_vs_Day7_WT"), each = 4),
#'   cell_type = rep(c("CD8 T cells", "B cells"), each = 2, times = 2),
#'   regulation = rep(c("Up", "Down"), times = 4),
#'   count = c(25, 10, 15, 5, 30, 12, 18, 8)
#' )
#'
#' # Create plot with triangles and numbers inside
#' gene_regulation_plot(
#'   data = sample_data,
#'   comparison_col = "comparison",
#'   cell_type_col = "cell_type",
#'   regulation_col = "regulation",
#'   count_col = "count",
#'   count_display = "inside",
#'   use_triangles = TRUE
#' )
#'
#' # Create plot with circles and numbers below
#' gene_regulation_plot(
#'   data = sample_data,
#'   comparison_col = "comparison",
#'   cell_type_col = "cell_type",
#'   regulation_col = "regulation",
#'   count_col = "count",
#'   count_display = "below",
#'   use_triangles = FALSE
#' )
#'
gene_regulation_plot <- function(data, 
                                comparison_col = "comparison_pattern", 
                                cell_type_col = "cell_type", 
                                regulation_col = "regulation", 
                                count_col = "count",
                                count_display = "none",
                                use_triangles = TRUE,
                                point_size_range = c(1, 10),
                                up_color = "deepskyblue",
                                down_color = "darkorange",
                                offset = 0.2,
                                title = "DEGs Across Cell Type and Comparisons",
                                x_angle = 45) {
  
  # Ensure required libraries are available
  require(ggplot2)
  require(dplyr)
  
  # Validate count_display parameter
  if (!count_display %in% c("none", "inside", "below", "both")) {
    stop("count_display must be one of: 'none', 'inside', 'below', or 'both'")
  }
  
  # Create a copy of the data to avoid modifying the original
  plot_data <- data
  
  # Rename columns to standardized names
  if (comparison_col != "comparison_pattern") {
    plot_data$comparison_pattern <- plot_data[[comparison_col]]
  }
  
  if (cell_type_col != "cell_type") {
    plot_data$cell_type <- plot_data[[cell_type_col]]
  }
  
  if (regulation_col != "regulation") {
    plot_data$regulation <- plot_data[[regulation_col]]
  }
  
  if (count_col != "count") {
    plot_data$count <- plot_data[[count_col]]
  }
  
  # Create numeric position variable for the x-axis
  plot_data <- plot_data %>%
    # Get unique comparisons
    mutate(comp_id = as.numeric(factor(comparison_pattern))) %>%
    # Add offset for regulation (Up/Down)
    mutate(x_pos = comp_id + ifelse(regulation == "Up", -offset, offset)) %>%
    # Convert cell_type to a factor
    mutate(cell_type = factor(cell_type))
  
  # Set shapes based on the use_triangles parameter
  if (use_triangles) {
    # Triangle shapes (up/down)
    shape_values <- c("Up" = 24, "Down" = 25)
    
    # Create the base plot with triangles
    p <- ggplot(plot_data, aes(
      x = x_pos, 
      y = cell_type, 
      size = count, 
      color = regulation,
      shape = regulation
    )) +
      geom_point(alpha = 0.9, stroke = 2, fill = "white") +
      scale_shape_manual(values = shape_values)
  } else {
    # Circle shape for both, color differentiates up/down
    # Create the base plot with circles
    p <- ggplot(plot_data, aes(
      x = x_pos, 
      y = cell_type, 
      size = count, 
      color = regulation
    )) +
      geom_point(alpha = 0.9, shape = 21, stroke = 2) +
      scale_color_manual(values = c("Up" = alpha(up_color, 0.3), "Down" = alpha(down_color, 0.3)))
  }
  
  # Add common aesthetic elements
  p <- p +
    scale_color_manual(values = c("Up" = up_color, "Down" = down_color)) +
    scale_size_continuous(range = point_size_range) +
    # Use unique comparison names for x-axis labels
    scale_x_continuous(
      breaks = unique(plot_data$comp_id),
      labels = unique(plot_data$comparison_pattern)
    ) +
    labs(
      title = title,
      x = "Comparison",
      y = "Cell Type",
      size = "Gene Count",
      color = "Regulation"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = x_angle, hjust = 1),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      legend.position = "right"
    )
  
  # Add shape to legend if using triangles
  if (use_triangles) {
    p <- p + labs(shape = "Regulation")
  } else {
    p <- p + labs(fill = "Regulation")
  }
  
  # Add text based on count_display parameter
  if (count_display %in% c("inside", "both")) {
    p <- p + geom_text(aes(label = count), 
                       size = 2.5, 
                       show.legend = FALSE,
                       color = "black")
  }
  
  if (count_display %in% c("below", "both")) {
    p <- p + geom_text(aes(label = count), 
                       nudge_y = -0.3,
                       size = 3, 
                       show.legend = FALSE)
  }
  
  return(p)
}