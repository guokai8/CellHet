# compareResponses.R

#' Compare Cell Type Responses to a Reference Cell Type
#'
#' @param deg_results Output from compareDEGs function
#' @param reference_cell_type Reference cell type to compare others against
#' @param plot_type Type of plot ("heatmap", "barplot", "radar")
#' @param similarity_metric Metric to calculate response similarity
#'
#' @return A plot showing how cell types respond differentially compared to reference
#' @export
compareResponses <- function(deg_results,
                             reference_cell_type,
                             plot_type = "heatmap",
                             similarity_metric = "jaccard") {

  # Extract DEG results
  deg_list <- deg_results$degs

  if (length(deg_list) == 0) {
    stop("No DEG results found in the input data")
  }

  # Extract cell types and comparisons from DEG results
  result_keys <- names(deg_list)
  result_parts <- strsplit(result_keys, "__")
  cell_types <- sapply(result_parts, function(x) x[1])
  comparisons <- sapply(result_parts, function(x) x[2])

  unique_cell_types <- unique(cell_types)
  unique_comparisons <- unique(comparisons)

  # Check if reference cell type exists
  if (!reference_cell_type %in% unique_cell_types) {
    stop(paste0("Reference cell type '", reference_cell_type, "' not found in results"))
  }

  # Calculate response similarity for each cell type compared to reference
  response_matrix <- matrix(0,
                            nrow = length(unique_cell_types),
                            ncol = length(unique_comparisons))
  rownames(response_matrix) <- unique_cell_types
  colnames(response_matrix) <- unique_comparisons

  # For each comparison, calculate similarity to reference
  for (comp in unique_comparisons) {
    # Get reference DEGs
    ref_key <- paste(reference_cell_type, comp, sep = "__")

    if (ref_key %in% names(deg_list)) {
      ref_genes <- deg_list[[ref_key]]$gene[deg_list[[ref_key]]$significant]

      # Compare to each other cell type
      for (ct in unique_cell_types) {
        ct_key <- paste(ct, comp, sep = "__")

        if (ct_key %in% names(deg_list)) {
          ct_genes <- deg_list[[ct_key]]$gene[deg_list[[ct_key]]$significant]

          # Calculate similarity
          similarity <- calculateFeatureSimilarity(ref_genes, ct_genes, similarity_metric)

          # Store in matrix (1 - similarity = divergence from reference)
          response_matrix[ct, comp] <- 1 - similarity
        }
      }
    }
  }

  # Create visualization based on plot_type
  if (plot_type == "heatmap") {
    # Create heatmap of response divergence
    if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      # Use ComplexHeatmap
      plot <- ComplexHeatmap::Heatmap(
        matrix = response_matrix,
        name = "Divergence",
        row_title = "Cell Type",
        column_title = "Comparison",
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 10),
        column_names_rot = 45,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid::grid.text(sprintf("%.2f", response_matrix[i, j]), x, y,
                          gp = grid::gpar(fontsize = 10))
        }
      )
    } else {
      # Fallback to geom_tile
      heatmap_data <- reshape2::melt(response_matrix,
                                     varnames = c("cell_type", "comparison"),
                                     value.name = "divergence")

      plot <- ggplot(heatmap_data, aes(x = comparison, y = cell_type, fill = divergence)) +
        geom_tile() +
        geom_text(aes(label = sprintf("%.2f", divergence))) +
        scale_fill_viridis_c(option = "magma") +
        theme_minimal() +
        labs(
          title = paste0("Cell Type Response Divergence from ", reference_cell_type),
          x = "Comparison",
          y = "Cell Type",
          fill = "Divergence"
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }

  } else if (plot_type == "barplot") {
    # Calculate average divergence for each cell type
    avg_divergence <- rowMeans(response_matrix, na.rm = TRUE)
    avg_data <- data.frame(
      cell_type = names(avg_divergence),
      divergence = avg_divergence
    )

    # Order by divergence
    avg_data <- avg_data[order(avg_data$divergence, decreasing = TRUE), ]

    # Create barplot
    plot <- ggplot(avg_data, aes(x = reorder(cell_type, -divergence), y = divergence)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      labs(
        title = paste0("Average Response Divergence from ", reference_cell_type),
        x = "Cell Type",
        y = "Divergence Score"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

  } else if (plot_type == "radar") {
    if (!requireNamespace("fmsb", quietly = TRUE)) {
      warning("Package 'fmsb' needed for radar plots. Falling back to heatmap.")
      return(compareResponses(deg_results, reference_cell_type, "heatmap", similarity_metric))
    }

    # Prepare data for radar plot
    radar_data <- t(response_matrix)

    # Add min and max rows required by fmsb
    radar_data <- rbind(rep(1, ncol(radar_data)), rep(0, ncol(radar_data)), radar_data)

    # Create radar plot
    par(mar = c(1, 1, 3, 1))
    plot <- fmsb::radarchart(
      radar_data,
      axistype = 1,
      pcol = rainbow(nrow(radar_data) - 2),
      plwd = 2,
      plty = 1,
      cglcol = "grey",
      cglty = 1,
      axislabcol = "grey",
      caxislabels = seq(0, 1, 0.25),
      cglwd = 0.8,
      title = paste0("Cell Type Response Divergence from ", reference_cell_type)
    )

    # Add legend
    legend("topright",
           legend = rownames(response_matrix),
           col = rainbow(nrow(radar_data) - 2),
           lty = 1,
           lwd = 2,
           cex = 0.8)

  } else {
    stop(paste0("Unsupported plot type: ", plot_type))
  }

  # Return results
  return(list(
    plot = plot,
    response_matrix = response_matrix,
    reference_cell_type = reference_cell_type
  ))
}
