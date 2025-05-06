#' Compare Cell Type Responses Across Conditions with Same Reference
#'
#' Visualizes how cell types respond when compared against the same reference condition.
#' Useful for identifying consistent patterns of cell type-specific responses.
#'
#' @param deg_results Output from compareDEGs function or findDifferentialGenes function
#' @param reference_group The reference group/condition used in all comparisons
#' @param cell_type_var Name of the column containing cell type annotations (default: "cell_type")
#' @param plot_type Type of visualization: "heatmap", "similarity", "trajectory" or "all"
#' @param min_genes Minimum number of DEGs required to include in visualization
#' @param cluster_cell_types Whether to cluster cell types by response similarity
#' @param cluster_conditions Whether to cluster conditions by response similarity
#' @param color_palette Color palette to use for visualization
#' @param show_numbers Whether to show count numbers on heatmap tiles
#' @param interactive Whether to create interactive plots (requires plotly)
#'
#' @return A list of plots or a combined plot object
#' @export
visualizeReferenceComparisons <- function(deg_results,
                                          reference_group,
                                          cell_type_var = "cell_type",
                                          plot_type = "heatmap",
                                          min_genes = 10,
                                          cluster_cell_types = TRUE,
                                          cluster_conditions = FALSE,
                                          color_palette = NULL,
                                          show_numbers = TRUE,
                                          interactive = FALSE) {

  # Check required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function. Please install it.")
  }

  # Extract data based on input type
  if (inherits(deg_results, "list")) {
    if ("all_degs" %in% names(deg_results)) {
      # Data from findDifferentialGenes
      all_degs <- deg_results$all_degs
    } else if ("degs" %in% names(deg_results)) {
      # Data from compareDEGs
      all_degs <- data.frame()

      for (result_key in names(deg_results$degs)) {
        # Extract cell type and comparison from key
        key_parts <- strsplit(result_key, "__")[[1]]
        cell_type <- key_parts[1]
        comparison <- key_parts[2]

        # Get DEG results
        degs <- deg_results$degs[[result_key]]

        # Add cell type and comparison columns
        degs$cell_type <- cell_type
        degs$comparison <- comparison

        # Split comparison into groups
        comp_parts <- strsplit(comparison, "_vs_")[[1]]
        degs$group1 <- comp_parts[1]
        degs$group2 <- comp_parts[2]

        # Append to combined data frame
        all_degs <- rbind(all_degs, degs)
      }
    } else {
      stop("Cannot find DEG results in the input object")
    }
  } else if (inherits(deg_results, "data.frame")) {
    # Input is already a data frame
    all_degs <- deg_results
  } else {
    stop("Input must be a data frame or a list containing DEG results")
  }

  # Filter for comparisons with the reference group
  all_degs <- all_degs[all_degs$group2 == reference_group | all_degs$group1 == reference_group, ]

  # Ensure consistent comparison orientation (reference is always group2)
  swap_indices <- which(all_degs$group1 == reference_group)
  if (length(swap_indices) > 0) {
    # Swap group labels
    tmp <- all_degs$group1[swap_indices]
    all_degs$group1[swap_indices] <- all_degs$group2[swap_indices]
    all_degs$group2[swap_indices] <- tmp

    # Invert log fold change
    all_degs$avg_log2FC[swap_indices] <- -all_degs$avg_log2FC[swap_indices]

    # Swap direction label
    all_degs$direction[swap_indices] <- ifelse(all_degs$direction[swap_indices] == "up", "down", "up")
  }

  # Update comparison label for consistency
  all_degs$comparison <- paste(all_degs$group1, "vs", all_degs$group2, sep = "_")

  # Check if we have data after filtering
  if (nrow(all_degs) == 0) {
    stop("No differential expression results found with the specified reference group")
  }

  # Get unique cell types and conditions
  unique_cell_types <- unique(all_degs$cell_type)
  unique_conditions <- unique(all_degs$group1)  # These are the non-reference conditions

  # Initialize result list
  result_plots <- list()

  # Function to calculate summary statistics for each cell type and condition
  calculate_summary_stats <- function() {
    # Create a matrix for:
    # 1. Number of up-regulated genes
    up_matrix <- matrix(0, nrow = length(unique_cell_types), ncol = length(unique_conditions))
    rownames(up_matrix) <- unique_cell_types
    colnames(up_matrix) <- unique_conditions

    # 2. Number of down-regulated genes
    down_matrix <- matrix(0, nrow = length(unique_cell_types), ncol = length(unique_conditions))
    rownames(down_matrix) <- unique_cell_types
    colnames(down_matrix) <- unique_conditions

    # 3. Average log fold change
    logfc_matrix <- matrix(0, nrow = length(unique_cell_types), ncol = length(unique_conditions))
    rownames(logfc_matrix) <- unique_cell_types
    colnames(logfc_matrix) <- unique_conditions

    # Populate matrices
    for (ct in unique_cell_types) {
      for (cond in unique_conditions) {
        # Filter for this cell type and condition
        subset_degs <- all_degs[all_degs$cell_type == ct & all_degs$group1 == cond & all_degs$significant == TRUE, ]

        if (nrow(subset_degs) > 0) {
          # Count up/down genes
          up_matrix[ct, cond] <- sum(subset_degs$direction == "up")
          down_matrix[ct, cond] <- sum(subset_degs$direction == "down")

          # Calculate average logFC
          logfc_matrix[ct, cond] <- mean(subset_degs$avg_log2FC)
        }
      }
    }

    # Calculate response ratio (up/down)
    response_ratio <- matrix(0, nrow = length(unique_cell_types), ncol = length(unique_conditions))
    rownames(response_ratio) <- unique_cell_types
    colnames(response_ratio) <- unique_conditions

    for (i in 1:nrow(response_ratio)) {
      for (j in 1:ncol(response_ratio)) {
        if (down_matrix[i, j] > 0) {
          response_ratio[i, j] <- up_matrix[i, j] / down_matrix[i, j]
        } else if (up_matrix[i, j] > 0) {
          response_ratio[i, j] <- Inf
        } else {
          response_ratio[i, j] <- NA
        }
      }
    }

    # Calculate total DEGs
    total_degs <- up_matrix + down_matrix

    # Calculate response score (combination of magnitude and direction)
    response_score <- matrix(0, nrow = length(unique_cell_types), ncol = length(unique_conditions))
    rownames(response_score) <- unique_cell_types
    colnames(response_score) <- unique_conditions

    for (i in 1:nrow(response_score)) {
      for (j in 1:ncol(response_score)) {
        if (total_degs[i, j] > 0) {
          # Calculate weighted score: sign(logFC) * log10(total_DEGs)
          response_score[i, j] <- sign(logfc_matrix[i, j]) * log10(total_degs[i, j] + 1)
        }
      }
    }

    # Return all matrices
    return(list(
      up = up_matrix,
      down = down_matrix,
      logfc = logfc_matrix,
      ratio = response_ratio,
      total = total_degs,
      response_score = response_score
    ))
  }

  # Calculate summary statistics
  summary_stats <- calculate_summary_stats()

  # Filter by minimum gene count
  for (ct in rownames(summary_stats$total)) {
    for (cond in colnames(summary_stats$total)) {
      if (summary_stats$total[ct, cond] < min_genes) {
        # Reset all values for this cell type and condition
        summary_stats$up[ct, cond] <- 0
        summary_stats$down[ct, cond] <- 0
        summary_stats$logfc[ct, cond] <- 0
        summary_stats$ratio[ct, cond] <- NA
        summary_stats$total[ct, cond] <- 0
        summary_stats$response_score[ct, cond] <- 0
      }
    }
  }

  # Create heatmap visualization
  if (plot_type %in% c("heatmap", "all")) {
    # Cluster cell types if requested
    if (cluster_cell_types) {
      # Use response_score matrix for clustering
      ct_dist <- stats::dist(summary_stats$response_score, method = "euclidean")
      ct_hclust <- stats::hclust(ct_dist, method = "ward.D2")
      ct_order <- rownames(summary_stats$response_score)[ct_hclust$order]
    } else {
      ct_order <- rownames(summary_stats$response_score)
    }

    # Cluster conditions if requested
    if (cluster_conditions) {
      # Transpose response_score matrix for clustering
      cond_dist <- stats::dist(t(summary_stats$response_score), method = "euclidean")
      cond_hclust <- stats::hclust(cond_dist, method = "ward.D2")
      cond_order <- colnames(summary_stats$response_score)[cond_hclust$order]
    } else {
      cond_order <- colnames(summary_stats$response_score)
    }

    # Create data frame for up-regulated genes
    up_df <- reshape2::melt(summary_stats$up,
                            varnames = c("cell_type", "condition"),
                            value.name = "up_genes")
    up_df$cell_type <- factor(up_df$cell_type, levels = ct_order)
    up_df$condition <- factor(up_df$condition, levels = cond_order)

    # Create data frame for down-regulated genes
    down_df <- reshape2::melt(summary_stats$down,
                              varnames = c("cell_type", "condition"),
                              value.name = "down_genes")

    # Create data frame for log fold change
    logfc_df <- reshape2::melt(summary_stats$logfc,
                               varnames = c("cell_type", "condition"),
                               value.name = "avg_logfc")

    # Create data frame for total DEGs
    total_df <- reshape2::melt(summary_stats$total,
                               varnames = c("cell_type", "condition"),
                               value.name = "total_degs")

    # Create data frame for response score
    score_df <- reshape2::melt(summary_stats$response_score,
                               varnames = c("cell_type", "condition"),
                               value.name = "response_score")

    # Combine data frames
    plot_df <- dplyr::left_join(up_df, down_df, by = c("cell_type", "condition"))
    plot_df <- dplyr::left_join(plot_df, logfc_df, by = c("cell_type", "condition"))
    plot_df <- dplyr::left_join(plot_df, total_df, by = c("cell_type", "condition"))
    plot_df <- dplyr::left_join(plot_df, score_df, by = c("cell_type", "condition"))

    # Calculate up/down ratio
    plot_df$ratio <- plot_df$up_genes / (plot_df$down_genes + 0.001)  # Add small value to avoid division by zero

    # Cap ratio for visualization
    plot_df$ratio_capped <- pmin(plot_df$ratio, 10)

    # Create heatmap of total DEGs
    total_heatmap <- ggplot2::ggplot(plot_df, ggplot2::aes(x = condition, y = cell_type, fill = total_degs)) +
      ggplot2::geom_tile(color = "white", size = 0.5) +
      ggplot2::scale_fill_viridis_c(name = "Total DEGs", option = "plasma") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Total DEGs vs", reference_group),
        x = "Condition",
        y = "Cell Type"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank()
      )

    if (show_numbers) {
      total_heatmap <- total_heatmap +
        ggplot2::geom_text(ggplot2::aes(label = total_degs), size = 3)
    }

    # Create heatmap of response score
    score_heatmap <- ggplot2::ggplot(plot_df, ggplot2::aes(x = condition, y = cell_type, fill = response_score)) +
      ggplot2::geom_tile(color = "white", size = 0.5) +
      ggplot2::scale_fill_gradient2(
        low = "blue", mid = "white", high = "red", midpoint = 0,
        name = "Response Score"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Cell Type Response Score vs", reference_group),
        subtitle = "Direction and magnitude of response (red = up, blue = down)",
        x = "Condition",
        y = "Cell Type"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank()
      )

    if (show_numbers) {
      score_heatmap <- score_heatmap +
        ggplot2::geom_text(ggplot2::aes(label = round(response_score, 2)), size = 3)
    }

    # Create heatmap of up/down ratio
    ratio_heatmap <- ggplot2::ggplot(plot_df, ggplot2::aes(x = condition, y = cell_type, fill = ratio_capped)) +
      ggplot2::geom_tile(color = "white", size = 0.5) +
      ggplot2::scale_fill_gradient(
        low = "lightblue", high = "darkred",
        name = "Up/Down Ratio",
        trans = "log10"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Up/Down Regulation Ratio vs", reference_group),
        x = "Condition",
        y = "Cell Type"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank()
      )

    if (show_numbers) {
      ratio_heatmap <- ratio_heatmap +
        ggplot2::geom_text(ggplot2::aes(label = round(ratio, 2)), size = 3)
    }

    # Make interactive if requested
    if (interactive && requireNamespace("plotly", quietly = TRUE)) {
      total_heatmap <- plotly::ggplotly(total_heatmap)
      score_heatmap <- plotly::ggplotly(score_heatmap)
      ratio_heatmap <- plotly::ggplotly(ratio_heatmap)
    }

    # Store heatmaps in result list
    result_plots$total_heatmap <- total_heatmap
    result_plots$score_heatmap <- score_heatmap
    result_plots$ratio_heatmap <- ratio_heatmap
  }

  # Create radar plot
  if (plot_type %in% c("radar", "all")) {
    if (!requireNamespace("ggradar", quietly = TRUE)) {
      warning("Package 'ggradar' required for radar plots. Skipping radar visualization.")
    } else {
      # Prepare data for radar plot
      # For each cell type, normalize response scores across conditions
      radar_data <- data.frame(cell_type = unique_cell_types)

      for (cond in unique_conditions) {
        radar_data[[cond]] <- scale(summary_stats$response_score[, cond])[, 1]
      }

      # Replace NAs with 0
      radar_data[is.na(radar_data)] <- 0

      # Scale to 0-1 range for radar plot
      for (i in 2:ncol(radar_data)) {
        col_min <- min(radar_data[, i])
        col_max <- max(radar_data[, i])

        if (col_max > col_min) {
          radar_data[, i] <- (radar_data[, i] - col_min) / (col_max - col_min)
        } else {
          radar_data[, i] <- 0.5  # Set to middle if all values are the same
        }
      }

      # Create radar plot for each cell type
      radar_plots <- list()

      for (ct in unique_cell_types) {
        ct_data <- radar_data[radar_data$cell_type == ct, ]

        # Skip if no data
        if (nrow(ct_data) == 0) next

        # Create radar plot
        radar_plot <- ggradar::ggradar(
          ct_data,
          values.radar = c("0", "0.5", "1"),
          grid.min = 0, grid.max = 1,
          grid.label.size = 3,
          axis.label.size = 3,
          legend.position = "none",
          plot.title = paste("Response Profile for", ct, "vs", reference_group)
        )

        radar_plots[[ct]] <- radar_plot
      }

      # Combine radar plots if possible
      if (requireNamespace("patchwork", quietly = TRUE) && length(radar_plots) > 0) {
        # Calculate appropriate layout
        n_plots <- length(radar_plots)
        n_cols <- min(3, n_plots)

        combined_radar <- patchwork::wrap_plots(
          radar_plots,
          ncol = n_cols
        ) +
          patchwork::plot_annotation(
            title = paste("Cell Type Response Profiles vs", reference_group),
            theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5))
          )

        result_plots$radar_plots <- combined_radar
      } else {
        # Store individual plots
        result_plots$radar_plots <- radar_plots
      }
    }
  }

  # Create cell type similarity visualization
  if (plot_type %in% c("similarity", "all")) {
    # Calculate similarity matrix between cell types based on response patterns
    similarity_matrix <- matrix(0,
                                nrow = length(unique_cell_types),
                                ncol = length(unique_cell_types))
    rownames(similarity_matrix) <- unique_cell_types
    colnames(similarity_matrix) <- unique_cell_types

    # Fill similarity matrix (using correlation of response scores)
    for (i in 1:length(unique_cell_types)) {
      for (j in i:length(unique_cell_types)) {
        ct1 <- unique_cell_types[i]
        ct2 <- unique_cell_types[j]

        if (i == j) {
          # Perfect similarity with self
          similarity_matrix[ct1, ct2] <- 1
        } else {
          # Calculate correlation between response patterns
          scores1 <- summary_stats$response_score[ct1, ]
          scores2 <- summary_stats$response_score[ct2, ]

          # Compute correlation, handle NA/NaN
          cor_val <- suppressWarnings(stats::cor(scores1, scores2, method = "pearson"))

          if (is.na(cor_val)) {
            similarity_matrix[ct1, ct2] <- similarity_matrix[ct2, ct1] <- 0
          } else {
            # Scale from -1:1 to 0:1
            similarity_matrix[ct1, ct2] <- similarity_matrix[ct2, ct1] <- (cor_val + 1) / 2
          }
        }
      }
    }

    # Create heatmap of cell type similarity
    similarity_df <- reshape2::melt(similarity_matrix,
                                    varnames = c("CellType1", "CellType2"),
                                    value.name = "Similarity")

    # Create heatmap
    similarity_heatmap <- ggplot2::ggplot(similarity_df,
                                          ggplot2::aes(x = CellType1, y = CellType2, fill = Similarity)) +
      ggplot2::geom_tile(color = "white", size = 0.5) +
      ggplot2::scale_fill_gradient(low = "white", high = "darkblue", limits = c(0, 1)) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Cell Type Response Similarity",
        subtitle = paste("Based on response patterns vs", reference_group),
        x = "Cell Type",
        y = "Cell Type"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank()
      )

    # Add similarity values
    similarity_heatmap <- similarity_heatmap +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Similarity)), size = 3)

    # Create network visualization if igraph is available
    if (requireNamespace("igraph", quietly = TRUE) &&
        requireNamespace("ggraph", quietly = TRUE)) {

      # Create network from similarity matrix
      # Only keep edges with similarity above threshold
      similarity_threshold <- 0.6  # Adjustable threshold

      # Create edge list
      edges <- data.frame()
      for (i in 1:nrow(similarity_matrix)) {
        for (j in i:ncol(similarity_matrix)) {
          if (i != j && similarity_matrix[i, j] >= similarity_threshold) {
            edges <- rbind(edges, data.frame(
              from = rownames(similarity_matrix)[i],
              to = colnames(similarity_matrix)[j],
              weight = similarity_matrix[i, j]
            ))
          }
        }
      }

      # Skip network if no edges
      if (nrow(edges) > 0) {
        # Create graph
        g <- igraph::graph_from_data_frame(edges, directed = FALSE,
                                           vertices = data.frame(name = unique_cell_types))

        # Calculate node size based on number of DEGs
        node_sizes <- rowSums(summary_stats$total)
        igraph::V(g)$size <- log10(node_sizes + 1) * 5

        # Create network plot
        network_plot <- ggraph::ggraph(g, layout = "fr") +
          ggraph::geom_edge_link(ggplot2::aes(width = weight, alpha = weight),
                                 edge_colour = "steelblue") +
          ggraph::geom_node_point(ggplot2::aes(size = size), color = "darkblue") +
          ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 4) +
          ggraph::scale_edge_width(range = c(1, 3)) +
          ggraph::scale_edge_alpha(range = c(0.6, 1)) +
          ggplot2::theme_void() +
          ggplot2::labs(
            title = "Cell Type Response Similarity Network",
            subtitle = paste("Based on response patterns vs", reference_group)
          ) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
            plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5),
            legend.position = "none"
          )

        result_plots$similarity_network <- network_plot
      }
    }

    result_plots$similarity_heatmap <- similarity_heatmap
  }

  # Create response trajectory visualization
  if (plot_type %in% c("trajectory", "all")) {
    # Calculate trajectory for each cell type across conditions
    # Use average log fold change as trajectory measure
    trajectory_data <- reshape2::melt(summary_stats$logfc,
                                      varnames = c("cell_type", "condition"),
                                      value.name = "avg_logfc")

    # Add up and down gene counts
    up_counts <- reshape2::melt(summary_stats$up,
                                varnames = c("cell_type", "condition"),
                                value.name = "up_genes")

    down_counts <- reshape2::melt(summary_stats$down,
                                  varnames = c("cell_type", "condition"),
                                  value.name = "down_genes")

    trajectory_data <- dplyr::left_join(trajectory_data, up_counts, by = c("cell_type", "condition"))
    trajectory_data <- dplyr::left_join(trajectory_data, down_counts, by = c("cell_type", "condition"))

    # Calculate total DEGs and set point size
    trajectory_data$total_degs <- trajectory_data$up_genes + trajectory_data$down_genes
    trajectory_data$point_size <- log10(trajectory_data$total_degs + 1) * 3

    # Create line plot of response trajectories
    trajectory_plot <- ggplot2::ggplot(trajectory_data,
                                       ggplot2::aes(x = condition, y = avg_logfc,
                                                    group = cell_type, color = cell_type)) +
      ggplot2::geom_line(size = 1) +
      ggplot2::geom_point(ggplot2::aes(size = point_size)) +
      ggplot2::scale_size_identity() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Response Trajectories vs", reference_group),
        x = "Condition",
        y = "Average Log2 Fold Change",
        color = "Cell Type"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "right"
      ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")

    # Make interactive if requested
    if (interactive && requireNamespace("plotly", quietly = TRUE)) {
      trajectory_plot <- plotly::ggplotly(trajectory_plot)
    }

    result_plots$trajectory_plot <- trajectory_plot

    # Create a faceted version with one panel per cell type
    facet_plot <- ggplot2::ggplot(trajectory_data,
                                  ggplot2::aes(x = condition, y = avg_logfc,
                                               group = 1, color = cell_type)) +
      ggplot2::geom_line(size = 1) +
      ggplot2::geom_point(ggplot2::aes(size = point_size)) +
      ggplot2::scale_size_identity() +
      ggplot2::facet_wrap(~ cell_type, scales = "free_y") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Cell Type-Specific Response Trajectories vs", reference_group),
        x = "Condition",
        y = "Average Log2 Fold Change"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.background = ggplot2::element_rect(fill = "lightgray"),
        strip.text = ggplot2::element_text(face = "bold")
      ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")

    result_plots$faceted_trajectory_plot <- facet_plot
  }

  # Return all plots or just the selected type
  if (plot_type == "all") {
    # Combine plots with patchwork if available
    if (requireNamespace("patchwork", quietly = TRUE) && !interactive) {
      main_plots <- list()

      # Add key plots to main display
      if ("score_heatmap" %in% names(result_plots)) {
        main_plots$score_heatmap <- result_plots$score_heatmap
      }

      if ("similarity_heatmap" %in% names(result_plots)) {
        main_plots$similarity_heatmap <- result_plots$similarity_heatmap
      }

      if ("trajectory_plot" %in% names(result_plots)) {
        main_plots$trajectory_plot <- result_plots$trajectory_plot
      }

      # Combine main plots
      if (length(main_plots) > 0) {
        combined_plot <- patchwork::wrap_plots(
          main_plots,
          ncol = 1
        ) +
          patchwork::plot_annotation(
            title = paste("Cell Type Response Patterns vs", reference_group),
            theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5))
          )

        result_plots$combined_plot <- combined_plot
      }
    }

    return(result_plots)
  } else {
    # Return only the requested plot type
    return(result_plots[grep(plot_type, names(result_plots))])
  }
}

#' Extract DEG data from different inputs
#' @keywords internal
extractDEGData <- function(deg_results) {
  if (inherits(deg_results, "list")) {
    if ("all_degs" %in% names(deg_results)) {
      return(deg_results$all_degs)
    } else if ("degs" %in% names(deg_results)) {
      all_degs <- data.frame()
      for (result_key in names(deg_results$degs)) {
        key_parts <- strsplit(result_key, "__")[[1]]
        cell_type <- key_parts[1]
        comparison <- key_parts[2]

        degs <- deg_results$degs[[result_key]]
        degs$cell_type <- cell_type
        degs$comparison <- comparison
        all_degs <- rbind(all_degs, degs)
      }
      return(all_degs)
    }
  } else if (inherits(deg_results, "data.frame")) {
    return(deg_results)
  }
  stop("Input must be a data frame or a list containing DEG results")
}

#' Filter DEGs by direction
#' @keywords internal
filterByDirection <- function(sig_degs, direction) {
  if (direction == "up") {
    if ("direction" %in% colnames(sig_degs)) {
      return(sig_degs[sig_degs$direction == "up", ])
    } else if ("avg_log2FC" %in% colnames(sig_degs)) {
      return(sig_degs[sig_degs$avg_log2FC > 0, ])
    }
  } else if (direction == "down") {
    if ("direction" %in% colnames(sig_degs)) {
      return(sig_degs[sig_degs$direction == "down", ])
    } else if ("avg_log2FC" %in% colnames(sig_degs)) {
      return(sig_degs[sig_degs$avg_log2FC < 0, ])
    }
  }
  return(sig_degs)  # Return all for "both" or if direction columns not found
}

#' Filter by cell types and comparisons
#' @keywords internal
filterByGroups <- function(sig_degs, cell_types, comparisons) {
  if (!is.null(cell_types)) {
    sig_degs <- sig_degs[sig_degs$cell_type %in% cell_types, ]
  }
  if (!is.null(comparisons)) {
    sig_degs <- sig_degs[sig_degs$comparison %in% comparisons, ]
  }
  return(sig_degs)
}

#' Create gene lists based on grouping
#' @keywords internal
createGeneLists <- function(sig_degs, gene_col, by_cell_type) {
  gene_lists <- list()

  if (by_cell_type) {
    # Group by cell type
    for (ct in unique(sig_degs$cell_type)) {
      ct_genes <- unique(sig_degs[sig_degs$cell_type == ct, gene_col])
      if (length(ct_genes) > 0) {
        gene_lists[[ct]] <- ct_genes
      }
    }
  } else {
    # Group by comparison
    for (comp in unique(sig_degs$comparison)) {
      comp_genes <- unique(sig_degs[sig_degs$comparison == comp, gene_col])
      if (length(comp_genes) > 0) {
        gene_lists[[comp]] <- comp_genes
      }
    }
  }

  return(gene_lists)
}
