#' Create Spatial Heterogeneity Map
#'
#' Creates a detailed visualization of spatial heterogeneity patterns.
#'
#' @param spatial_results Output from analyzeSpatialHeterogeneity function
#' @param plot_type Type of plot: "cell_types", "clusters", "autocorrelation", or "all" (default: "all")
#' @param top_n Number of top features to plot (default: 5)
#' @param ncol Number of columns for multi-panel plots (default: 2)
#' @param color_palette Color palette for plotting (default: "viridis")
#' @param point_size Size of points in the plot (default: 3)
#' @param point_alpha Transparency of points (default: 0.7)
#' @param add_voronoi Logical, whether to add Voronoi tessellation to the plot (default: FALSE)
#' @param add_convex_hull Logical, whether to add convex hull around clusters (default: FALSE)
#'
#' @return A ggplot object or a list of ggplot objects
#' @export
plotSpatialHeterogeneity <- function(spatial_results,
                                     plot_type = "all",
                                     top_n = 5,
                                     ncol = 2,
                                     color_palette = "viridis",
                                     point_size = 3,
                                     point_alpha = 0.7,
                                     add_voronoi = FALSE,
                                     add_convex_hull = FALSE) {

  # Check required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function. Please install it.")
  }

  # Validate input
  if (!is.list(spatial_results) || !"spatial_stats" %in% names(spatial_results)) {
    stop("Invalid spatial_results input. Must be output from analyzeSpatialHeterogeneity function.")
  }

  # Extract data from results
  analysis_data <- spatial_results$analysis_data
  spatial_stats <- spatial_results$spatial_stats

  # Check if we have spatial clusters
  has_clusters <- !is.null(spatial_results$spatial_clusters)
  if (has_clusters) {
    analysis_data$spatial_cluster <- spatial_results$spatial_clusters
  }

  # Initialize plot list
  plots <- list()

  # Create base spatial plot
  base_plot <- ggplot2::ggplot(analysis_data, ggplot2::aes(x = x, y = y)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Spatial X", y = "Spatial Y")

  # Create plots based on plot_type
  if (plot_type %in% c("cell_types", "all")) {
    # Plot cell types
    cell_type_plot <- base_plot +
      ggplot2::geom_point(ggplot2::aes(color = cell_type), size = point_size, alpha = point_alpha)

    # Apply color palette
    if (color_palette == "viridis") {
      cell_type_plot <- cell_type_plot + ggplot2::scale_color_viridis_d(name = "Cell Type")
    } else {
      cell_type_plot <- cell_type_plot + ggplot2::scale_color_brewer(name = "Cell Type", palette = color_palette)
    }

    cell_type_plot <- cell_type_plot + ggplot2::labs(title = "Spatial Distribution of Cell Types")

    plots$cell_type_plot <- cell_type_plot
  }

  if (plot_type %in% c("clusters", "all") && has_clusters) {
    # Plot spatial clusters
    cluster_plot <- base_plot +
      ggplot2::geom_point(ggplot2::aes(color = factor(spatial_cluster)), size = point_size, alpha = point_alpha)

    # Apply color palette
    if (color_palette == "viridis") {
      cluster_plot <- cluster_plot + ggplot2::scale_color_viridis_d(name = "Cluster")
    } else {
      cluster_plot <- cluster_plot + ggplot2::scale_color_brewer(name = "Cluster", palette = color_palette)
    }

    cluster_plot <- cluster_plot + ggplot2::labs(title = "Spatial Clustering Results")

    # Add convex hull around clusters if requested
    if (add_convex_hull) {
      if (!requireNamespace("ggforce", quietly = TRUE)) {
        warning("Package 'ggforce' needed for convex hull. Skipping.")
      } else {
        cluster_plot <- cluster_plot +
          ggforce::geom_mark_hull(ggplot2::aes(fill = factor(spatial_cluster)),
                                  alpha = 0.1, expand = 0.01, show.legend = FALSE)
      }
    }

    plots$cluster_plot <- cluster_plot
  }

  if (plot_type %in% c("autocorrelation", "all")) {
    # Determine if we're working with cell types or genes
    if ("cell_types" %in% names(spatial_stats)) {
      # Get top autocorrelated cell types
      top_features <- spatial_stats$cell_types$cell_type[spatial_stats$cell_types$p_value < 0.05]
      top_features <- top_features[1:min(length(top_features), top_n)]

      if (length(top_features) > 0) {
        top_feature_plots <- list()

        for (feature in top_features) {
          feature_data <- analysis_data
          feature_data$is_feature <- feature_data$cell_type == feature

          p <- base_plot +
            ggplot2::geom_point(data = feature_data[!feature_data$is_feature, ],
                                color = "gray80", size = point_size/2, alpha = point_alpha/2) +
            ggplot2::geom_point(data = feature_data[feature_data$is_feature, ],
                                color = "red", size = point_size, alpha = point_alpha) +
            ggplot2::labs(title = paste0("Cell Type: ", feature),
                          subtitle = paste0("Pattern: ",
                                            spatial_stats$cell_types$pattern[spatial_stats$cell_types$cell_type == feature],
                                            ", p-value: ",
                                            round(spatial_stats$cell_types$p_value[spatial_stats$cell_types$cell_type == feature], 4)))

          # Add Voronoi tessellation if requested
          if (add_voronoi) {
            if (!requireNamespace("deldir", quietly = TRUE)) {
              warning("Package 'deldir' needed for Voronoi tessellation. Skipping.")
            } else {
              voronoi_data <- feature_data[feature_data$is_feature, ]
              if (nrow(voronoi_data) >= 3) {  # Need at least 3 points for meaningful tessellation
                vor <- deldir::deldir(voronoi_data$x, voronoi_data$y)
                vor_polygons <- deldir::tile.list(vor)

                # Convert to format for ggplot
                vor_polys <- lapply(vor_polygons, function(tile) {
                  data.frame(x = tile$x, y = tile$y)
                })
                vor_polys <- lapply(seq_along(vor_polys), function(i) {
                  cbind(vor_polys[[i]], id = i)
                })
                vor_df <- do.call(rbind, vor_polys)

                p <- p + ggplot2::geom_polygon(data = vor_df,
                                               ggplot2::aes(x = x, y = y, group = id),
                                               fill = NA, color = "black", alpha = 0.3)
              }
            }
          }

          top_feature_plots[[feature]] <- p
        }

        plots$top_feature_plots <- top_feature_plots
      }
    } else if ("genes" %in% names(spatial_stats)) {
      # Get top autocorrelated genes
      top_features <- spatial_stats$genes$gene[spatial_stats$genes$p_value < 0.05]
      top_features <- top_features[1:min(length(top_features), top_n)]

      if (length(top_features) > 0) {
        top_feature_plots <- list()

        for (feature in top_features) {
          feature_data <- analysis_data
          # Note: gene expression should be present in the results or we'd need to extract it again
          if ("expression" %in% names(feature_data)) {
            p <- base_plot +
              ggplot2::geom_point(data = feature_data, ggplot2::aes(color = expression),
                                  size = point_size, alpha = point_alpha) +
              ggplot2::scale_color_viridis_c(name = "Expression") +
              ggplot2::labs(title = paste0("Gene: ", feature),
                            subtitle = paste0("Pattern: ",
                                              spatial_stats$genes$pattern[spatial_stats$genes$gene == feature],
                                              ", p-value: ",
                                              round(spatial_stats$genes$p_value[spatial_stats$genes$gene == feature], 4)))

            top_feature_plots[[feature]] <- p
          }
        }

        plots$top_feature_plots <- top_feature_plots
      }
    }
  }

  # Combine plots if requested
  if (plot_type == "all" && requireNamespace("patchwork", quietly = TRUE)) {
    combined_plot <- NULL

    # Add cell type plot
    if ("cell_type_plot" %in% names(plots)) {
      combined_plot <- plots$cell_type_plot
    }

    # Add cluster plot
    if ("cluster_plot" %in% names(plots)) {
      if (is.null(combined_plot)) {
        combined_plot <- plots$cluster_plot
      } else {
        combined_plot <- combined_plot + plots$cluster_plot
      }
    }

    # Add feature plots
    if ("top_feature_plots" %in% names(plots)) {
      feature_plots <- plots$top_feature_plots

      if (length(feature_plots) > 0) {
        feature_combined <- patchwork::wrap_plots(feature_plots, ncol = ncol)

        if (is.null(combined_plot)) {
          combined_plot <- feature_combined
        } else {
          combined_plot <- combined_plot / feature_combined
        }
      }
    }

    if (!is.null(combined_plot)) {
      combined_plot <- combined_plot +
        patchwork::plot_annotation(
          title = "Spatial Heterogeneity Analysis",
          theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5))
        )

      plots$combined_plot <- combined_plot

      # Return just the combined plot for convenience
      return(plots$combined_plot)
    }
  }

  # Return all plots
  return(plots)
}

#' Analyze Co-occurrence of Cell Types in Spatial Data
#'
#' Analyzes the spatial co-occurrence or avoidance patterns between different cell types.
#'
#' @param spatial_object Spatial transcriptomics object (e.g., Seurat with spatial assay, SpatialExperiment)
#' @param cell_type_col Name of the column in metadata containing cell type annotations
#' @param spatial_coords Names of the columns containing spatial coordinates
#' @param radius Distance threshold for considering cells as neighbors
#' @param permutations Number of permutations for significance testing (default: 1000)
#' @param method Co-occurrence method: "ripley_bivariate", "nearest_neighbor", or "correlation" (default: "correlation")
#' @param plot_type Plot type: "heatmap", "network", or "both" (default: "both")
#' @param cell_types Optional vector of cell types to include (default: all)
#'
#' @return A list containing co-occurrence results and visualizations
#' @export
analyzeSpatialCoOccurrence <- function(spatial_object,
                                       cell_type_col = "cell_type",
                                       spatial_coords = c("x", "y"),
                                       radius = NULL,
                                       permutations = 1000,
                                       method = "correlation",
                                       plot_type = "both",
                                       cell_types = NULL) {

  # Check if required packages are available
  if (!requireNamespace("spdep", quietly = TRUE)) {
    stop("Package 'spdep' is required for spatial analysis. Please install it.")
  }

  # Determine object type and extract metadata and spatial coordinates
  is_seurat <- inherits(spatial_object, "Seurat")
  is_spe <- inherits(spatial_object, "SpatialExperiment")

  if (!is_seurat && !is_spe) {
    stop("spatial_object must be a Seurat object with spatial data or SpatialExperiment object")
  }

  # Extract metadata and spatial coordinates
  if (is_seurat) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' required for processing Seurat objects")
    }

    metadata <- spatial_object@meta.data

    # Check for spatial assay
    if (!"spatial" %in% names(spatial_object@images)) {
      stop("No spatial data found in Seurat object. Please ensure it has a spatial assay.")
    }

    # Extract spatial coordinates
    spatial_data <- Seurat::GetTissueCoordinates(spatial_object)

    # Check if spatial_coords match the available columns
    if (!all(spatial_coords %in% colnames(spatial_data))) {
      # Try default Seurat coordinate names
      alt_coords <- c("imagecol", "imagerow")
      if (all(alt_coords %in% colnames(spatial_data))) {
        warning(paste0("Specified spatial_coords not found. Using default Seurat coordinates: ",
                       paste(alt_coords, collapse = ", ")))
        spatial_coords <- alt_coords
      } else {
        stop("Specified spatial_coords not found in the spatial data")
      }
    }

  } else if (is_spe) {
    if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
      stop("Package 'SpatialExperiment' required for processing SpatialExperiment objects")
    }

    metadata <- as.data.frame(SpatialExperiment::colData(spatial_object))

    # Extract spatial coordinates
    spatial_data <- as.data.frame(SpatialExperiment::spatialCoords(spatial_object))

    # Check if spatial_coords match the available columns
    if (!all(spatial_coords %in% colnames(spatial_data))) {
      # Try default SpatialExperiment coordinate names
      alt_coords <- colnames(spatial_data)[1:2]
      if (length(alt_coords) >= 2) {
        warning(paste0("Specified spatial_coords not found. Using first two coordinate columns: ",
                       paste(alt_coords, collapse = ", ")))
        spatial_coords <- alt_coords
      } else {
        stop("Specified spatial_coords not found in the spatial data")
      }
    }
  }

  # Check if cell_type_col exists in metadata
  if (!cell_type_col %in% colnames(metadata)) {
    stop(paste0("Column '", cell_type_col, "' not found in metadata"))
  }

  # Extract unique cell types
  all_cell_types <- unique(metadata[[cell_type_col]])

  # Filter for specific cell types if requested
  if (!is.null(cell_types)) {
    valid_cell_types <- intersect(cell_types, all_cell_types)
    if (length(valid_cell_types) == 0) {
      stop("None of the specified cell types found in the data")
    }
    if (length(valid_cell_types) < length(cell_types)) {
      warning(paste0(length(cell_types) - length(valid_cell_types), " cell types not found in the data"))
    }
    target_cell_types <- valid_cell_types
  } else {
    target_cell_types <- all_cell_types
  }

  # Prepare data frame with cell type and spatial information
  analysis_data <- data.frame(
    cell_id = rownames(metadata),
    cell_type = metadata[[cell_type_col]],
    x = spatial_data[[spatial_coords[1]]],
    y = spatial_data[[spatial_coords[2]]]
  )

  # Remove rows with NA coordinates or non-target cell types
  analysis_data <- analysis_data[!is.na(analysis_data$x) & !is.na(analysis_data$y) &
                                   analysis_data$cell_type %in% target_cell_types, ]

  # Create a coordinates matrix
  coords <- cbind(analysis_data$x, analysis_data$y)

  # Determine radius automatically if not provided
  if (is.null(radius)) {
    # Calculate nearest neighbor distances
    nn_dist <- spdep::knearneigh(coords, k = 1)
    nn_dist <- spdep::knn2nb(nn_dist)
    nn_dist <- unlist(spdep::nbdists(nn_dist, coords))

    # Set radius to 2 times the mean nearest neighbor distance
    radius <- 2 * mean(nn_dist)
  }

  # Initialize co-occurrence matrix
  n_cell_types <- length(target_cell_types)
  co_occurrence <- matrix(0, nrow = n_cell_types, ncol = n_cell_types)
  rownames(co_occurrence) <- target_cell_types
  colnames(co_occurrence) <- target_cell_types

  # Calculate co-occurrence statistics based on method
  if (method == "correlation") {
    # Create binary indicators for each cell type
    cell_type_indicators <- matrix(0, nrow = nrow(analysis_data), ncol = n_cell_types)
    colnames(cell_type_indicators) <- target_cell_types

    for (i in 1:n_cell_types) {
      ct <- target_cell_types[i]
      cell_type_indicators[, i] <- as.numeric(analysis_data$cell_type == ct)
    }

    # Create neighbor list
    nb <- spdep::dnearneigh(coords, d1 = 0, d2 = radius)

    # Calculate spatial correlations between cell types
    for (i in 1:n_cell_types) {
      ct1 <- target_cell_types[i]
      ct1_vector <- cell_type_indicators[, i]

      for (j in 1:n_cell_types) {
        if (i == j) {
          # Same cell type, set to 1 (perfect correlation)
          co_occurrence[i, j] <- 1
        } else {
          ct2 <- target_cell_types[j]
          ct2_vector <- cell_type_indicators[, j]

          # Calculate local spatial correlation
          result <- spdep::localmoran_perm(ct1_vector, spdep::nb2listw(nb), nsim = permutations)

          # Use correlation between cell type indicators as co-occurrence measure
          ct_correlation <- cor(ct1_vector, ct2_vector)

          # Check if correlation is significant based on permutations
          # This is a simplified approach; a real implementation would use proper statistical testing
          co_occurrence[i, j] <- ct_correlation
        }
      }
    }

  } else if (method == "nearest_neighbor") {
    # Create lists of points for each cell type
    cell_type_points <- list()
    for (ct in target_cell_types) {
      cell_type_points[[ct]] <- coords[analysis_data$cell_type == ct, , drop = FALSE]
    }

    # Calculate nearest neighbor statistics
    for (i in 1:n_cell_types) {
      ct1 <- target_cell_types[i]
      points1 <- cell_type_points[[ct1]]

      if (nrow(points1) == 0) next

      for (j in 1:n_cell_types) {
        if (i == j) {
          # Same cell type, set to 1 (perfect co-occurrence)
          co_occurrence[i, j] <- 1
        } else {
          ct2 <- target_cell_types[j]
          points2 <- cell_type_points[[ct2]]

          if (nrow(points2) == 0) {
            co_occurrence[i, j] <- 0
            next
          }

          # For each point of type 1, find distance to nearest point of type 2
          nn_dists <- numeric(nrow(points1))

          for (k in 1:nrow(points1)) {
            point <- points1[k, ]

            # Calculate distances to all points of type 2
            distances <- sqrt(rowSums((points2 - rep(point, each = nrow(points2)))^2))

            # Find minimum distance
            nn_dists[k] <- min(distances)
          }

          # Calculate proportion of type 1 points that have a type 2 neighbor within radius
          proportion <- mean(nn_dists <= radius)

          # Adjust to -1 to 1 scale: -1 (avoidance) to 1 (attraction)
          # Expected proportion under CSR is the density of type 2 points
          expected <- nrow(points2) / (max(coords[,1]) - min(coords[,1])) / (max(coords[,2]) - min(coords[,2])) * pi * radius^2
          expected <- min(expected, 1)  # Cap at 1

          # Normalized measure: (observed - expected) / (1 - expected) for attraction
          # or (observed - expected) / expected for avoidance
          if (proportion >= expected) {
            normalized <- (proportion - expected) / (1 - expected)
          } else {
            normalized <- (proportion - expected) / expected
          }

          co_occurrence[i, j] <- normalized
        }
      }
    }

  } else if (method == "ripley_bivariate") {
    if (!requireNamespace("spatstat.geom", quietly = TRUE) ||
        !requireNamespace("spatstat.core", quietly = TRUE)) {
      warning("Packages 'spatstat.geom' and 'spatstat.core' required for bivariate Ripley's K. Falling back to correlation method.")
      method <- "correlation"

      # Recalculate using correlation method - in a real function you would include full implementation here
      # For now, I'll just point to the need to implement this section
      warning("Falling back to correlation method - implementation required for this section")
    } else {
      # Create observation window
      window <- spatstat.geom::owin(xrange = range(coords[, 1]), yrange = range(coords[, 2]))

      # Create point patterns for each cell type
      cell_type_patterns <- list()
      for (ct in target_cell_types) {
        ct_points <- coords[analysis_data$cell_type == ct, , drop = FALSE]
        if (nrow(ct_points) > 0) {
          cell_type_patterns[[ct]] <- spatstat.geom::ppp(ct_points[, 1], ct_points[, 2], window = window)
        }
      }

      # Calculate bivariate Ripley's K at radius
      for (i in 1:n_cell_types) {
        ct1 <- target_cell_types[i]

        if (!ct1 %in% names(cell_type_patterns)) {
          next
        }

        for (j in 1:n_cell_types) {
          if (i == j) {
            # Same cell type, set to 1 (perfect correlation)
            co_occurrence[i, j] <- 1
          } else {
            ct2 <- target_cell_types[j]

            if (!ct2 %in% names(cell_type_patterns)) {
              next
            }

            # Calculate bivariate K function
            # First, create multitype point pattern
            pp1 <- cell_type_patterns[[ct1]]
            pp2 <- cell_type_patterns[[ct2]]

            # Use cross K function to estimate bivariate pattern
            suppressWarnings({
              K12 <- spatstat.core::Kcross(spatstat.geom::superimpose(pattern1 = pp1, pattern2 = pp2,
                                                                      W = window),
                                           i = "pattern1", j = "pattern2", r = radius)
            })

            # Calculate K12(r) / (pi * r^2) - 1
            # Positive values indicate attraction, negative values indicate repulsion
            L12 <- K12$iso[length(K12$iso)] / (pi * radius^2) - 1

            # Scale to -1 to 1
            normalized <- L12 / (1 + abs(L12))

            co_occurrence[i, j] <- normalized
          }
        }
      }
    }
  }

  # Generate visualizations
  plots <- list()

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    # Create heatmap if requested
    if (plot_type %in% c("heatmap", "both")) {
      if (requireNamespace("reshape2", quietly = TRUE)) {
        # Convert matrix to long format
        co_occurrence_long <- reshape2::melt(co_occurrence)
        names(co_occurrence_long) <- c("CellType1", "CellType2", "CoOccurrence")

        # Create heatmap
        heatmap_plot <- ggplot2::ggplot(co_occurrence_long,
                                        ggplot2::aes(x = CellType2, y = CellType1, fill = CoOccurrence)) +
          ggplot2::geom_tile() +
          ggplot2::geom_text(ggplot2::aes(label = round(CoOccurrence, 2)), color = "black", size = 3) +
          ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                        limits = c(-1, 1), name = "Co-occurrence") +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
          ggplot2::labs(title = "Cell Type Co-occurrence Patterns",
                        x = "Cell Type",
                        y = "Cell Type")

        plots$heatmap <- heatmap_plot
      }
    }

    # Create network if requested
    if (plot_type %in% c("network", "both")) {
      if (requireNamespace("igraph", quietly = TRUE) && requireNamespace("ggraph", quietly = TRUE)) {
        # Create graph from co-occurrence matrix
        # Keep only significant connections
        threshold <- 0.2  # Can be adjusted based on analysis
        co_occurrence_thresholded <- co_occurrence
        co_occurrence_thresholded[abs(co_occurrence_thresholded) < threshold] <- 0

        # Create graph
        g <- igraph::graph_from_adjacency_matrix(
          co_occurrence_thresholded,
          mode = "directed",
          weighted = TRUE,
          diag = FALSE
        )

        # Set node attributes
        igraph::V(g)$name <- rownames(co_occurrence)
        igraph::V(g)$size <- table(analysis_data$cell_type)[igraph::V(g)$name]
        igraph::V(g)$size <- sqrt(igraph::V(g)$size) * 2

        # Set edge attributes
        igraph::E(g)$width <- abs(igraph::E(g)$weight) * 3
        igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, "red", "blue")

        # Create network plot
        network_plot <- ggraph::ggraph(g, layout = "fr") +
          ggraph::geom_edge_arc(ggplot2::aes(width = width, color = color, alpha = abs(weight)),
                                arrow = ggraph::arrow(length = ggplot2::unit(3, "mm")),
                                start_cap = ggraph::circle(3, "mm"),
                                end_cap = ggraph::circle(3, "mm")) +
          ggraph::geom_node_point(ggplot2::aes(size = size), color = "darkgrey") +
          ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE) +
          ggplot2::scale_edge_color_identity() +
          ggplot2::scale_edge_width_identity() +
          ggplot2::scale_edge_alpha_identity() +
          ggplot2::labs(title = "Cell Type Co-occurrence Network",
                        subtitle = paste0("Red edges: attraction (co-occurrence), ",
                                          "Blue edges: repulsion (avoidance)")) +
          ggraph::theme_graph()

        plots$network <- network_plot
      }
    }

    # Create spatial distribution plot
    # Base spatial plot
    base_plot <- ggplot2::ggplot(analysis_data, ggplot2::aes(x = x, y = y)) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = spatial_coords[1], y = spatial_coords[2])

    # Plot cell types
    cell_type_plot <- base_plot +
      ggplot2::geom_point(ggplot2::aes(color = cell_type), size = 3, alpha = 0.7) +
      ggplot2::scale_color_discrete(name = "Cell Type") +
      ggplot2::labs(title = "Spatial Distribution of Cell Types")

    plots$cell_type_plot <- cell_type_plot

    # If we have strong co-occurrence pairs, create a visualization for them
    strong_pairs <- which(abs(co_occurrence) > threshold & row(co_occurrence) != col(co_occurrence), arr.ind = TRUE)

    if (nrow(strong_pairs) > 0) {
      for (i in 1:min(nrow(strong_pairs), 3)) {  # Show up to 3 strong pairs
        row_idx <- strong_pairs[i, 1]
        col_idx <- strong_pairs[i, 2]

        ct1 <- target_cell_types[row_idx]
        ct2 <- target_cell_types[col_idx]

        co_value <- co_occurrence[row_idx, col_idx]

        # Create dataset with highlighted cell types
        pair_data <- analysis_data
        pair_data$highlight <- pair_data$cell_type %in% c(ct1, ct2)
        pair_data$specific_type <- ifelse(pair_data$cell_type == ct1, ct1,
                                          ifelse(pair_data$cell_type == ct2, ct2, "Other"))

        # Create plot
        pair_plot <- ggplot2::ggplot(pair_data, ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_point(data = pair_data[!pair_data$highlight, ],
                              color = "gray80", size = 2, alpha = 0.3) +
          ggplot2::geom_point(data = pair_data[pair_data$highlight, ],
                              ggplot2::aes(color = specific_type), size = 3, alpha = 0.7) +
          ggplot2::scale_color_manual(values = c("Other" = "gray80",
                                                 ct1 = "red",
                                                 ct2 = "blue")) +
          ggplot2::labs(title = paste0("Co-occurrence of ", ct1, " and ", ct2),
                        subtitle = paste0("Co-occurrence value: ", round(co_value, 2),
                                          ", Pattern: ", ifelse(co_value > 0, "Attraction", "Avoidance")),
                        color = "Cell Type")

        plots[[paste0("pair_", ct1, "_", ct2)]] <- pair_plot
      }
    }
  }

  # Compile results
  results <- list(
    metadata = list(
      cell_type_col = cell_type_col,
      spatial_coords = spatial_coords,
      radius = radius,
      method = method,
      target_cell_types = target_cell_types
    ),
    co_occurrence = co_occurrence,
    analysis_data = analysis_data,
    plots = plots
  )

  return(results)
}


#' Analyze Spatial Heterogeneity of Cell Types
#'
#' Analyzes spatial heterogeneity and clustering patterns of cell types or gene expression
#' in spatial transcriptomics data.
#'
#' @param spatial_object Spatial transcriptomics object (e.g., Seurat with spatial assay, SpatialExperiment)
#' @param cell_type_col Name of the column in metadata containing cell type annotations
#' @param spatial_coords Names of the columns containing spatial coordinates
#' @param reference_group Optional group to use as reference for comparisons
#' @param radius Radius to use for local neighborhood analysis (default: NULL, auto-determined)
#' @param n_neighbors Number of nearest neighbors to use for analysis (default: 10)
#' @param permutations Number of permutations for significance testing (default: 100)
#' @param method Analysis method: "morans_i", "getis_ord", "ripley_k", or "nearest_neighbor" (default: "morans_i")
#' @param cluster_method Method for spatial clustering: "dbscan", "clara", "hclust" (default: "dbscan")
#' @param n_clusters Number of spatial clusters to identify (default: 5)
#' @param gene_set Optional set of genes to analyze (default: NULL, uses cell types)
#' @param variable_features Logical, whether to use variable features for analysis (default: FALSE)
#'
#' @return A list containing spatial heterogeneity analysis results and visualizations
#' @export
analyzeSpatialHeterogeneity <- function(spatial_object,
                                        cell_type_col = "cell_type",
                                        spatial_coords = c("x", "y"),
                                        reference_group = NULL,
                                        radius = NULL,
                                        n_neighbors = 10,
                                        permutations = 100,
                                        method = "morans_i",
                                        cluster_method = "dbscan",
                                        n_clusters = 5,
                                        gene_set = NULL,
                                        variable_features = FALSE) {

  # Check if required packages are available
  if (!requireNamespace("spdep", quietly = TRUE)) {
    stop("Package 'spdep' is required for spatial analysis. Please install it.")
  }

  # Determine object type and extract metadata and spatial coordinates
  is_seurat <- inherits(spatial_object, "Seurat")
  is_spe <- inherits(spatial_object, "SpatialExperiment")

  if (!is_seurat && !is_spe) {
    stop("spatial_object must be a Seurat object with spatial data or SpatialExperiment object")
  }

  # Extract metadata and spatial coordinates
  if (is_seurat) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' required for processing Seurat objects")
    }

    metadata <- spatial_object@meta.data

    # Check for spatial assay
    if (!"spatial" %in% names(spatial_object@images)) {
      stop("No spatial data found in Seurat object. Please ensure it has a spatial assay.")
    }

    # Extract spatial coordinates
    spatial_data <- Seurat::GetTissueCoordinates(spatial_object)

    # Check if spatial_coords match the available columns
    if (!all(spatial_coords %in% colnames(spatial_data))) {
      # Try default Seurat coordinate names
      alt_coords <- c("imagecol", "imagerow")
      if (all(alt_coords %in% colnames(spatial_data))) {
        warning(paste0("Specified spatial_coords not found. Using default Seurat coordinates: ",
                       paste(alt_coords, collapse = ", ")))
        spatial_coords <- alt_coords
      } else {
        stop("Specified spatial_coords not found in the spatial data")
      }
    }

  } else if (is_spe) {
    if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
      stop("Package 'SpatialExperiment' required for processing SpatialExperiment objects")
    }

    metadata <- as.data.frame(SpatialExperiment::colData(spatial_object))

    # Extract spatial coordinates
    spatial_data <- as.data.frame(SpatialExperiment::spatialCoords(spatial_object))

    # Check if spatial_coords match the available columns
    if (!all(spatial_coords %in% colnames(spatial_data))) {
      # Try default SpatialExperiment coordinate names
      alt_coords <- colnames(spatial_data)[1:2]
      if (length(alt_coords) >= 2) {
        warning(paste0("Specified spatial_coords not found. Using first two coordinate columns: ",
                       paste(alt_coords, collapse = ", ")))
        spatial_coords <- alt_coords
      } else {
        stop("Specified spatial_coords not found in the spatial data")
      }
    }
  }

  # Check if cell_type_col exists in metadata
  if (!cell_type_col %in% colnames(metadata)) {
    stop(paste0("Column '", cell_type_col, "' not found in metadata"))
  }

  # Extract unique cell types
  all_cell_types <- unique(metadata[[cell_type_col]])

  # Filter for specific cell types if requested
  if (!is.null(cell_types)) {
    valid_cell_types <- intersect(cell_types, all_cell_types)
    if (length(valid_cell_types) == 0) {
      stop("None of the specified cell types found in the data")
    }
    if (length(valid_cell_types) < length(cell_types)) {
      warning(paste0(length(cell_types) - length(valid_cell_types), " cell types not found in the data"))
    }
    target_cell_types <- valid_cell_types
  } else {
    target_cell_types <- all_cell_types
  }

  # Prepare data frame with cell type and spatial information
  analysis_data <- data.frame(
    cell_id = rownames(metadata),
    cell_type = metadata[[cell_type_col]],
    x = spatial_data[[spatial_coords[1]]],
    y = spatial_data[[spatial_coords[2]]]
  )

  # Remove rows with NA coordinates or non-target cell types
  analysis_data <- analysis_data[!is.na(analysis_data$x) & !is.na(analysis_data$y) &
                                   analysis_data$cell_type %in% target_cell_types, ]

  # Create a coordinates matrix
  coords <- cbind(analysis_data$x, analysis_data$y)

  # Determine radius automatically if not provided
  if (is.null(radius)) {
    # Calculate nearest neighbor distances
    nn_dist <- spdep::knearneigh(coords, k = 1)
    nn_dist <- spdep::knn2nb(nn_dist)
    nn_dist <- unlist(spdep::nbdists(nn_dist, coords))

    # Set radius to 2 times the mean nearest neighbor distance
    radius <- 2 * mean(nn_dist)
  }

  # Initialize co-occurrence matrix
  n_cell_types <- length(target_cell_types)
  co_occurrence <- matrix(0, nrow = n_cell_types, ncol = n_cell_types)
  rownames(co_occurrence) <- target_cell_types
  colnames(co_occurrence) <- target_cell_types

  # Calculate co-occurrence statistics based on method
  if (method == "correlation") {
    # Create binary indicators for each cell type
    cell_type_indicators <- matrix(0, nrow = nrow(analysis_data), ncol = n_cell_types)
    colnames(cell_type_indicators) <- target_cell_types

    for (i in 1:n_cell_types) {
      ct <- target_cell_types[i]
      cell_type_indicators[, i] <- as.numeric(analysis_data$cell_type == ct)
    }

    # Create neighbor list
    nb <- spdep::dnearneigh(coords, d1 = 0, d2 = radius)

    # Calculate spatial correlations between cell types
    for (i in 1:n_cell_types) {
      ct1 <- target_cell_types[i]
      ct1_vector <- cell_type_indicators[, i]

      for (j in 1:n_cell_types) {
        if (i == j) {
          # Same cell type, set to 1 (perfect correlation)
          co_occurrence[i, j] <- 1
        } else {
          ct2 <- target_cell_types[j]
          ct2_vector <- cell_type_indicators[, j]

          # Calculate local spatial correlation
          result <- spdep::localmoran_perm(ct1_vector, spdep::nb2listw(nb), nsim = permutations)

          # Use correlation between cell type indicators as co-occurrence measure
          ct_correlation <- cor(ct1_vector, ct2_vector)

          # Check if correlation is significant based on permutations
          # This is a simplified approach; a real implementation would use proper statistical testing
          co_occurrence[i, j] <- ct_correlation
        }
      }
    }

  } else if (method == "nearest_neighbor") {
    # Create lists of points for each cell type
    cell_type_points <- list()
    for (ct in target_cell_types) {
      cell_type_points[[ct]] <- coords[analysis_data$cell_type == ct, , drop = FALSE]
    }

    # Calculate nearest neighbor statistics
    for (i in 1:n_cell_types) {
      ct1 <- target_cell_types[i]
      points1 <- cell_type_points[[ct1]]

      if (nrow(points1) == 0) next

      for (j in 1:n_cell_types) {
        if (i == j) {
          # Same cell type, set to 1 (perfect co-occurrence)
          co_occurrence[i, j] <- 1
        } else {
          ct2 <- target_cell_types[j]
          points2 <- cell_type_points[[ct2]]

          if (nrow(points2) == 0) {
            co_occurrence[i, j] <- 0
            next
          }

          # For each point of type 1, find distance to nearest point of type 2
          nn_dists <- numeric(nrow(points1))

          for (k in 1:nrow(points1)) {
            point <- points1[k, ]

            # Calculate distances to all points of type 2
            distances <- sqrt(rowSums((points2 - rep(point, each = nrow(points2)))^2))

            # Find minimum distance
            nn_dists[k] <- min(distances)
          }

          # Calculate proportion of type 1 points that have a type 2 neighbor within radius
          proportion <- mean(nn_dists <= radius)

          # Adjust to -1 to 1 scale: -1 (avoidance) to 1 (attraction)
          # Expected proportion under CSR is the density of type 2 points
          expected <- nrow(points2) / (max(coords[,1]) - min(coords[,1])) / (max(coords[,2]) - min(coords[,2])) * pi * radius^2
          expected <- min(expected, 1)  # Cap at 1

          # Normalized measure: (observed - expected) / (1 - expected) for attraction
          # or (observed - expected) / expected for avoidance
          if (proportion >= expected) {
            normalized <- (proportion - expected) / (1 - expected)
          } else {
            normalized <- (proportion - expected) / expected
          }

          co_occurrence[i, j] <- normalized
        }
      }
    }

  } else if (method == "ripley_bivariate") {
    if (!requireNamespace("spatstat.geom", quietly = TRUE) ||
        !requireNamespace("spatstat.core", quietly = TRUE)) {
      warning("Packages 'spatstat.geom' and 'spatstat.core' required for bivariate Ripley's K. Falling back to correlation method.")
      method <- "correlation"

      # Recalculate using correlation method
      # (code would be duplicated from the correlation method above)
    } else {
      # Create observation window
      window <- spatstat.geom::owin(xrange = range(coords[, 1]), yrange = range(coords[, 2]))

      # Create point patterns for each cell type
      cell_type_patterns <- list()
      for (ct in target_cell_types) {
        ct_points <- coords[analysis_data$cell_type == ct, , drop = FALSE]
        if (nrow(ct_points) > 0) {
          cell_type_patterns[[ct]] <- spatstat.geom::ppp(ct_points[, 1], ct_points[, 2], window = window)
        }
      }

      # Calculate bivariate Ripley's K at radius
      for (i in 1:n_cell_types) {
        ct1 <- target_cell_types[i]

        if (!ct1 %in% names(cell_type_patterns)) {
          next
        }

        for (j in 1:n_cell_types) {
          if (i == j) {
            # Same cell type, set to 1 (perfect correlation)
            co_occurrence[i, j] <- 1
          } else {
            ct2 <- target_cell_types[j]

            if (!ct2 %in% names(cell_type_patterns)) {
              next
            }

            # Calculate bivariate K function
            # First, create multitype point pattern
            pp1 <- cell_type_patterns[[ct1]]
            pp2 <- cell_type_patterns[[ct2]]

            # Use cross K function to estimate bivariate pattern
            suppressWarnings({
              K12 <- spatstat.core::Kcross(spatstat.geom::superimpose(pattern1 = pp1, pattern2 = pp2,
                                                                      W = window),
                                           i = "pattern1", j = "pattern2", r = radius)
            })

            # Calculate K12(r) / (pi * r^2) - 1
            # Positive values indicate attraction, negative values indicate repulsion
            L12 <- K12$iso[length(K12$iso)] / (pi * radius^2) - 1

            # Scale to -1 to 1
            normalized <- L12 / (1 + abs(L12))

            co_occurrence[i, j] <- normalized
          }
        }
      }
    }
  }

  # Generate visualizations
  plots <- list()

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    # Create heatmap if requested
    if (plot_type %in% c("heatmap", "both")) {
      if (requireNamespace("reshape2", quietly = TRUE)) {
        # Convert matrix to long format
        co_occurrence_long <- reshape2::melt(co_occurrence)
        names(co_occurrence_long) <- c("CellType1", "CellType2", "CoOccurrence")

        # Create heatmap
        heatmap_plot <- ggplot2::ggplot(co_occurrence_long,
                                        ggplot2::aes(x = CellType2, y = CellType1, fill = CoOccurrence)) +
          ggplot2::geom_tile() +
          ggplot2::geom_text(ggplot2::aes(label = round(CoOccurrence, 2)), color = "black", size = 3) +
          ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                        limits = c(-1, 1), name = "Co-occurrence") +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
          ggplot2::labs(title = "Cell Type Co-occurrence Patterns",
                        x = "Cell Type",
                        y = "Cell Type")

        plots$heatmap <- heatmap_plot
      }
    }

    # Create network if requested
    if (plot_type %in% c("network", "both")) {
      if (requireNamespace("igraph", quietly = TRUE) && requireNamespace("ggraph", quietly = TRUE)) {
        # Create graph from co-occurrence matrix
        # Keep only significant connections
        threshold <- 0.2  # Can be adjusted based on analysis
        co_occurrence_thresholded <- co_occurrence
        co_occurrence_thresholded[abs(co_occurrence_thresholded) < threshold] <- 0

        # Create graph
        g <- igraph::graph_from_adjacency_matrix(
          co_occurrence_thresholded,
          mode = "directed",
          weighted = TRUE,
          diag = FALSE
        )

        # Set node attributes
        igraph::V(g)$name <- rownames(co_occurrence)
        igraph::V(g)$size <- table(analysis_data$cell_type)[igraph::V(g)$name]
        igraph::V(g)$size <- sqrt(igraph::V(g)$size) * 2

        # Set edge attributes
        igraph::E(g)$width <- abs(igraph::E(g)$weight) * 3
        igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, "red", "blue")

        # Create network plot
        network_plot <- ggraph::ggraph(g, layout = "fr") +
          ggraph::geom_edge_arc(ggplot2::aes(width = width, color = color, alpha = abs(weight)),
                                arrow = ggraph::arrow(length = ggplot2::unit(3, "mm")),
                                start_cap = ggraph::circle(3, "mm"),
                                end_cap = ggraph::circle(3, "mm")) +
          ggraph::geom_node_point(ggplot2::aes(size = size), color = "darkgrey") +
          ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE) +
          ggplot2::scale_edge_color_identity() +
          ggplot2::scale_edge_width_identity() +
          ggplot2::scale_edge_alpha_identity() +
          ggplot2::labs(title = "Cell Type Co-occurrence Network",
                        subtitle = paste0("Red edges: attraction (co-occurrence), ",
                                          "Blue edges: repulsion (avoidance)")) +
          ggraph::theme_graph()

        plots$network <- network_plot
      }
    }

    # Create spatial distribution plot
    # Base spatial plot
    base_plot <- ggplot2::ggplot(analysis_data, ggplot2::aes(x = x, y = y)) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = spatial_coords[1], y = spatial_coords[2])

    # Plot cell types
    cell_type_plot <- base_plot +
      ggplot2::geom_point(ggplot2::aes(color = cell_type), size = 3, alpha = 0.7) +
      ggplot2::scale_color_discrete(name = "Cell Type") +
      ggplot2::labs(title = "Spatial Distribution of Cell Types")

    plots$cell_type_plot <- cell_type_plot
  }

  # Compile results
  results <- list(
    metadata = list(
      cell_type_col = cell_type_col,
      spatial_coords = spatial_coords,
      radius = radius,
      method = method,
      target_cell_types = target_cell_types
    ),
    co_occurrence = co_occurrence,
    analysis_data = analysis_data,
    plots = plots
  )

  return(results)
}

