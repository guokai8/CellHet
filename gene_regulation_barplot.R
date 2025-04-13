#' Create a Barplot of Gene Regulation by Cell Type
#'
#' Creates a barplot showing up and down regulated genes for each cell type
#' across different comparisons.
#'
#' @param data A data frame containing DEG results
#' @param comparison_col Column name for the comparison categories
#' @param cell_type_col Column name for the cell types
#' @param up_col Column name for up-regulated gene counts
#' @param down_col Column name for down-regulated gene counts
#' @param up_color Color for up-regulated genes (default = "darkorange")
#' @param down_color Color for down-regulated genes (default = "deepskyblue")
#' @param show_counts Logical. Whether to show count numbers on bars (default = TRUE)
#' @param count_size Size of count numbers (default = 3)
#' @param count_color Color of count numbers (default = "black")
#' @param horizontal Logical. If TRUE, creates horizontal bars (default = FALSE)
#' @param bar_width Width of the bars (default = 0.7)
#' @param facet_by Parameter to facet the plot by: "comparison", "cell_type", or "none" (default = "comparison")
#' @param facet_scales Scale parameter for facets: "free", "free_x", "free_y", or "fixed" (default = "free_y")
#' @param facet_ncol Number of columns in the facet grid (default = NULL)
#' @param stacked Logical. If TRUE, creates stacked bars (default = FALSE)
#' @param title Plot title (default = "Up and Down Regulated Genes by Cell Type")
#' @param title_size Size of the title (default = 12)
#' @param x_lab Label for x-axis (default = "Cell Type" or "Gene Count" if horizontal)
#' @param y_lab Label for y-axis (default = "Gene Count" or "Cell Type" if horizontal)
#' @param legend_title Title for the legend (default = "Regulation")
#' @param x_angle Angle for x-axis text (default = 45)
#'
#' @return A ggplot object
#'
gene_regulation_barplot <- function(data,
                                    comparison_col = "comparison_pattern",
                                    cell_type_col = "cell_type",
                                    up_col = "up_genes",
                                    down_col = "down_genes",
                                    up_color = "darkorange",
                                    down_color = "deepskyblue",
                                    show_counts = TRUE,
                                    count_size = 3,
                                    count_color = "black",
                                    horizontal = FALSE,
                                    bar_width = 0.7,
                                    facet_by = "comparison",
                                    facet_scales = "free_y",
                                    facet_ncol = NULL,
                                    stacked = FALSE,
                                    title = "Up and Down Regulated Genes by Cell Type",
                                    title_size = 12,
                                    x_lab = NULL,
                                    y_lab = NULL,
                                    legend_title = "Regulation",
                                    x_angle = 45) {
  
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  
  # Create a copy of the data to avoid modifying the original
  plot_data <- data
  
  # Validate column names
  col_names <- colnames(plot_data)
  if (!comparison_col %in% col_names) stop(paste("Column", comparison_col, "not found in data"))
  if (!cell_type_col %in% col_names) stop(paste("Column", cell_type_col, "not found in data"))
  if (!up_col %in% col_names) stop(paste("Column", up_col, "not found in data"))
  if (!down_col %in% col_names) stop(paste("Column", down_col, "not found in data"))
  
  # Validate facet_by parameter
  if (!facet_by %in% c("comparison", "cell_type", "none")) {
    stop("facet_by must be one of: 'comparison', 'cell_type', or 'none'")
  }
  
  # Reshape data to long format for plotting
  plot_data_long <- plot_data %>%
    pivot_longer(
      cols = c(up_col, down_col),
      names_to = "regulation",
      values_to = "count"
    ) %>%
    mutate(
      regulation = ifelse(regulation == up_col, "Up", "Down"),
      # Make sure cell_type and comparison_pattern are factors to preserve order
      !!sym(cell_type_col) := factor(!!sym(cell_type_col)),
      !!sym(comparison_col) := factor(!!sym(comparison_col))
    )
  
  # Default axis labels based on horizontal orientation
  if (is.null(x_lab)) {
    x_lab <- if (horizontal) "Gene Count" else "Cell Type"
  }
  if (is.null(y_lab)) {
    y_lab <- if (horizontal) "Cell Type" else "Gene Count"
  }
  
  # Determine position type (dodge or stack)
  position_type <- if (stacked) "stack" else "dodge"
  
  # Create the base plot
  if (horizontal) {
    # Horizontal bars
    p <- ggplot(plot_data_long) +
      geom_col(
        aes(
          y = !!sym(cell_type_col),
          x = count,
          fill = regulation
        ),
        position = position_type,
        width = bar_width
      )
  } else {
    # Vertical bars
    p <- ggplot(plot_data_long) +
      geom_col(
        aes(
          x = !!sym(cell_type_col),
          y = count,
          fill = regulation
        ),
        position = position_type,
        width = bar_width
      )
  }
  
  # Add counts to the bars if requested
  if (show_counts) {
    if (horizontal) {
      # For horizontal bars
      if (stacked) {
        # For stacked bars, position text in the middle of each segment
        p <- p + geom_text(
          aes(
            y = !!sym(cell_type_col),
            x = count/2,
            label = count,
            group = regulation
          ),
          position = position_stack(),
          size = count_size,
          color = count_color
        )
      } else {
        # For dodged bars, position text at the end of each bar
        p <- p + geom_text(
          aes(
            y = !!sym(cell_type_col),
            x = count,
            label = count,
            group = regulation
          ),
          position = position_dodge(width = bar_width),
          hjust = -0.2,
          size = count_size,
          color = count_color
        )
      }
    } else {
      # For vertical bars
      if (stacked) {
        # For stacked bars, position text in the middle of each segment
        p <- p + geom_text(
          aes(
            x = !!sym(cell_type_col),
            y = count/2,
            label = count,
            group = regulation
          ),
          position = position_stack(),
          size = count_size,
          color = count_color
        )
      } else {
        # For dodged bars, position text at the top of each bar
        p <- p + geom_text(
          aes(
            x = !!sym(cell_type_col),
            y = count,
            label = count,
            group = regulation
          ),
          position = position_dodge(width = bar_width),
          vjust = -0.2,
          size = count_size,
          color = count_color
        )
      }
    }
  }
  
  # Add fill colors and legend
  p <- p + scale_fill_manual(
    values = c("Up" = up_color, "Down" = down_color),
    name = legend_title
  )
  
  # Add facets if requested
  if (facet_by != "none") {
    if (facet_by == "comparison") {
      p <- p + facet_wrap(
        vars(!!sym(comparison_col)),
        scales = facet_scales,
        ncol = facet_ncol
      )
    } else if (facet_by == "cell_type") {
      p <- p + facet_wrap(
        vars(!!sym(cell_type_col)),
        scales = facet_scales,
        ncol = facet_ncol
      )
    }
  }
  
  # Complete the plot styling
  p <- p + labs(
    title = title,
    x = x_lab,
    y = y_lab
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = title_size, face = "bold"),
    axis.text.x = if (!horizontal) element_text(angle = x_angle, hjust = 1) else element_text(),
    legend.position = "right",
    panel.grid.major.x = if (horizontal) element_line(color = "grey90") else element_blank(),
    panel.grid.major.y = if (!horizontal) element_line(color = "grey90") else element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold")
  )
  
  return(p)
}

# Example usage:
# 1. Basic vertical barplot faceted by comparison
# gene_regulation_barplot(
#   deg_results,
#   comparison_col = "comparison_pattern",
#   cell_type_col = "cell_type",
#   up_col = "up_genes",
#   down_col = "down_genes",
#   facet_by = "comparison"
# )

# 2. Horizontal barplot for better readability with many cell types
# gene_regulation_barplot(
#   deg_results,
#   horizontal = TRUE,
#   facet_by = "comparison",
#   facet_ncol = 2
# )

# 3. Stacked barplot to show total regulated genes
# gene_regulation_barplot(
#   deg_results,
#   stacked = TRUE,
#   facet_by = "comparison"
# )

# 4. Barplot grouped by cell types to compare across conditions
# gene_regulation_barplot(
#   deg_results,
#   facet_by = "cell_type",
#   facet_scales = "free",
#   facet_ncol = 3
# )