#' Summarize DEG Comparison Results
#'
#' @param deg_results Output from compareDEGs function
#' @param plot_type Type of summary plot ("barplot", "heatmap", "dotplot")
#' @param direction Filter for gene direction ("up", "down", "both")
#' @param interactive Make plot interactive (requires plotly)
#'
#' @return A ggplot object showing the number of DEGs per comparison
#' @export
summarizeDEGs <- function(deg_results,
                          plot_type = "barplot",
                          direction = "both",
                          interactive = FALSE) {

  # Extract summary data
  summary_data <- deg_results$summary

  if (nrow(summary_data) == 0) {
    stop("No comparison results available in the input data")
  }

  # Filter by direction if needed
  if (direction == "up") {
    summary_data$n_genes <- summary_data$n_up
  } else if (direction == "down") {
    summary_data$n_genes <- summary_data$n_down
  } else {
    summary_data$n_genes <- summary_data$n_significant
  }

  # Create appropriate visualization
  if (plot_type == "barplot") {
    # Create barplot
    p <- ggplot2::ggplot(summary_data, ggplot2::aes(x = comparison, y = n_genes, fill = cell_type)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Number of",
                      ifelse(direction == "both", "significant", direction),
                      "DEGs by comparison"),
        x = "Comparison",
        y = "Number of DEGs",
        fill = "Cell Type"
      ) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  } else if (plot_type == "heatmap") {
    # Create heatmap
    heatmap_data <- reshape2::dcast(summary_data,
                                    cell_type ~ comparison,
                                    value.var = "n_genes")

    # Convert to matrix
    rownames(heatmap_data) <- heatmap_data$cell_type
    heatmap_data$cell_type <- NULL
    heatmap_matrix <- as.matrix(heatmap_data)

    # Create heatmap
    if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      # Use ComplexHeatmap
      p <- ComplexHeatmap::Heatmap(
        matrix = heatmap_matrix,
        name = "DEGs",
        row_title = "Cell Type",
        column_title = "Comparison",
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 10),
        column_names_rot = 45,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid::grid.text(sprintf("%d", heatmap_matrix[i, j]), x, y,
                          gp = grid::gpar(fontsize = 10))
        }
      )
    } else {
      # Fallback to geom_tile
      heatmap_data_long <- reshape2::melt(heatmap_matrix,
                                          varnames = c("cell_type", "comparison"),
                                          value.name = "n_genes")

      p <- ggplot2::ggplot(heatmap_data_long, ggplot2::aes(x = comparison, y = cell_type, fill = n_genes)) +
        ggplot2::geom_tile() +
        ggplot2::geom_text(ggplot2::aes(label = n_genes)) +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste("Number of",
                        ifelse(direction == "both", "significant", direction),
                        "DEGs by comparison"),
          x = "Comparison",
          y = "Cell Type",
          fill = "DEGs"
        ) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    }

  } else if (plot_type == "dotplot") {
    # Create dotplot
    p <- ggplot2::ggplot(summary_data,
                         ggplot2::aes(x = comparison, y = cell_type, size = n_genes, color = cell_type)) +
      ggplot2::geom_point() +
      ggplot2::scale_size_continuous(range = c(1, 10)) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Number of",
                      ifelse(direction == "both", "significant", direction),
                      "DEGs by comparison"),
        x = "Comparison",
        y = "Cell Type",
        size = "DEGs",
        color = "Cell Type"
      ) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  } else {
    stop(paste0("Unsupported plot type: ", plot_type))
  }

  # Make interactive if requested
  if (interactive) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("Package 'plotly' needed for interactive plots. Returning static plot instead.")
    } else if (plot_type != "heatmap" || !requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      p <- plotly::ggplotly(p)
    } else {
      warning("Interactive mode not supported for ComplexHeatmap. Returning static plot.")
    }
  }

  return(p)
}
