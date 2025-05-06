#' Create Comprehensive Visualizations for Multi-Modal Integration
#'
#' Generates a comprehensive set of visualizations for exploring the results of
#' multi-modal integration, including dimensionality reduction plots, feature expression
#' patterns, cluster comparisons, and quality metrics.
#'
#' @param integration_results Output from integrateMultiModalData function
#' @param color_by Vector of metadata columns to use for coloring points (default: "modality", "cell_type")
#' @param features Vector of features to visualize expression patterns
#' @param n_features Number of top variable features to visualize if features not specified
#' @param plot_types Vector of plot types to generate (default: all available plots)
#' @param split_by Optional metadata column to split visualization by
#' @param dim_reduction Dimensionality reduction to use for visualization (default: from integration)
#' @param interactive Logical, whether to create interactive plots using plotly (default: FALSE)
#' @param ncol Number of columns for arranging multiple plots (default: 2)
#' @param point_size Size of points in scatterplots (default: 1)
#' @param point_alpha Transparency of points (default: 0.7)
#' @param return_plots Logical, whether to return plot objects (default: TRUE)
#'
#' @return A list of plot objects and optionally an integrated plot panel
#' @export
visualizeIntegration <- function(integration_results,
                                 color_by = c("modality", "cell_type"),
                                 features = NULL,
                                 n_features = 5,
                                 plot_types = c("dimplot", "feature", "cluster", "metrics", "density"),
                                 split_by = NULL,
                                 dim_reduction = NULL,
                                 interactive = FALSE,
                                 ncol = 2,
                                 point_size = 1,
                                 point_alpha = 0.7,
                                 return_plots = TRUE) {

  # Extract integrated data and embeddings
  integrated_data <- integration_results$integration_results
  embeddings <- integration_results$integration_embeddings

  # Check if input is valid
  if (is.null(integrated_data) || is.null(embeddings)) {
    stop("No integration results or embeddings found")
  }

  # Determine object type
  is_seurat <- inherits(integrated_data, "Seurat")
  is_sce <- inherits(integrated_data, "SingleCellExperiment")

  if (!is_seurat && !is_sce) {
    stop("Integrated data must be a Seurat or SingleCellExperiment object")
  }

  # Extract metadata
  if (is_seurat) {
    metadata <- integrated_data@meta.data
  } else {
    metadata <- as.data.frame(SingleCellExperiment::colData(integrated_data))
  }

  # Validate color_by columns
  invalid_cols <- setdiff(color_by, colnames(metadata))
  if (length(invalid_cols) > 0) {
    warning(paste0("Columns not found in metadata: ", paste(invalid_cols, collapse = ", ")))
    color_by <- intersect(color_by, colnames(metadata))
  }

  if (length(color_by) == 0) {
    stop("No valid columns found for coloring")
  }

  # Validate split_by
  if (!is.null(split_by) && !split_by %in% colnames(metadata)) {
    warning(paste0("Split column '", split_by, "' not found in metadata. Ignoring."))
    split_by <- NULL
  }

  # Set dimensionality reduction
  if (is.null(dim_reduction)) {
    if (is_seurat) {
      # Use first available reduction
      reductions <- names(integrated_data@reductions)
      if ("umap" %in% reductions) {
        dim_reduction <- "umap"
      } else if ("tsne" %in% reductions) {
        dim_reduction <- "tsne"
      } else {
        dim_reduction <- reductions[1]
      }
    } else {
      # For SCE
      reductions <- SingleCellExperiment::reducedDimNames(integrated_data)
      if ("UMAP" %in% reductions) {
        dim_reduction <- "UMAP"
      } else if ("TSNE" %in% reductions) {
        dim_reduction <- "TSNE"
      } else {
        dim_reduction <- reductions[1]
      }
    }
  }

  # Extract embeddings for visualization
  if (is_seurat) {
    if (dim_reduction %in% names(integrated_data@reductions)) {
      plot_embeddings <- Seurat::Embeddings(integrated_data, reduction = dim_reduction)
    } else {
      warning(paste0("Reduction '", dim_reduction, "' not found. Using integration embeddings."))
      plot_embeddings <- embeddings
    }
  } else {
    if (dim_reduction %in% SingleCellExperiment::reducedDimNames(integrated_data)) {
      plot_embeddings <- SingleCellExperiment::reducedDim(integrated_data, dim_reduction)
    } else {
      warning(paste0("Reduction '", dim_reduction, "' not found. Using integration embeddings."))
      plot_embeddings <- embeddings
    }
  }

  # Prepare data frame for plotting
  plot_df <- data.frame(
    row.names = rownames(plot_embeddings),
    Dim1 = plot_embeddings[, 1],
    Dim2 = plot_embeddings[, 2]
  )

  # Add metadata columns
  for (col in colnames(metadata)) {
    plot_df[[col]] <- metadata[rownames(plot_df), col]
  }

  # Initialize results list
  plot_results <- list(
    plots = list(),
    combined_plot = NULL
  )

  # Check for required packages
  req_pkg <- c("ggplot2")
  if (interactive) {
    req_pkg <- c(req_pkg, "plotly")
  }

  missing_pkg <- req_pkg[!sapply(req_pkg, requireNamespace, quietly = TRUE)]
  if (length(missing_pkg) > 0) {
    warning(paste0("Missing required packages: ", paste(missing_pkg, collapse = ", ")))
    if ("ggplot2" %in% missing_pkg) {
      return(NULL)  # Can't continue without ggplot2
    }
    if ("plotly" %in% missing_pkg) {
      interactive <- FALSE  # Fall back to static plots
    }
  }

  # Generate dimension reduction plots colored by metadata
  if ("dimplot" %in% plot_types) {
    dimplots <- list()

    for (color_col in color_by) {
      # Create base plot
      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Dim1, y = Dim2, color = .data[[color_col]])) +
        ggplot2::geom_point(size = point_size, alpha = point_alpha) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste0(dim_reduction, " Plot Colored by ", color_col),
          x = paste0(dim_reduction, "1"),
          y = paste0(dim_reduction, "2"),
          color = color_col
        )

      # Add split if requested
      if (!is.null(split_by)) {
        p <- p + ggplot2::facet_wrap(~ .data[[split_by]])
      }

      # Make interactive if requested
      if (interactive) {
        p <- plotly::ggplotly(p)
      }

      # Store plot
      dimplots[[color_col]] <- p
    }

    plot_results$plots$dimplots <- dimplots
  }

  # Generate feature expression plots
  if ("feature" %in% plot_types) {
    # Determine features to plot
    plot_features <- features

    if (is.null(plot_features)) {
      # Find variable features if not specified
      if (is_seurat) {
        var_features <- Seurat::VariableFeatures(integrated_data)
        if (length(var_features) > 0) {
          plot_features <- var_features[1:min(n_features, length(var_features))]
        } else {
          # Get assay data and compute variance
          expr_data <- as.matrix(Seurat::GetAssayData(integrated_data, slot = "data"))
          gene_vars <- apply(expr_data, 1, stats::var)
          plot_features <- names(sort(gene_vars, decreasing = TRUE)[1:min(n_features, length(gene_vars))])
        }
      } else {
        # For SCE
        if ("modelGeneVar" %in% names(SummarizedExperiment::rowData(integrated_data))) {
          # Use existing variance modeling
          var_genes <- rownames(integrated_data)[order(SummarizedExperiment::rowData(integrated_data)$modelGeneVar$bio, decreasing = TRUE)]
          plot_features <- var_genes[1:min(n_features, length(var_genes))]
        } else {
          # Compute variance
          if ("logcounts" %in% SummarizedExperiment::assayNames(integrated_data)) {
            expr_data <- SummarizedExperiment::assay(integrated_data, "logcounts")
          } else {
            expr_data <- SummarizedExperiment::assay(integrated_data, 1)
          }
          gene_vars <- apply(expr_data, 1, stats::var)
          plot_features <- names(sort(gene_vars, decreasing = TRUE)[1:min(n_features, length(gene_vars))])
        }
      }
    }

    # Extract expression data for selected features
    if (is_seurat) {
      expr_data <- as.matrix(Seurat::GetAssayData(integrated_data, slot = "data"))
    } else {
      if ("logcounts" %in% SummarizedExperiment::assayNames(integrated_data)) {
        expr_data <- SummarizedExperiment::assay(integrated_data, "logcounts")
      } else {
        expr_data <- SummarizedExperiment::assay(integrated_data, 1)
      }
    }

    # Filter features to those present in the data
    plot_features <- intersect(plot_features, rownames(expr_data))

    if (length(plot_features) == 0) {
      warning("No features found for visualization")
    } else {
      # Create feature plots
      feature_plots <- list()

      for (feature in plot_features) {
        # Add expression to plot data
        plot_df[[feature]] <- expr_data[feature, rownames(plot_df)]

        # Create plot
        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Dim1, y = Dim2, color = .data[[feature]])) +
          ggplot2::geom_point(size = point_size, alpha = point_alpha) +
          ggplot2::scale_color_viridis_c() +
          ggplot2::theme_minimal() +
          ggplot2::labs(
            title = paste0("Expression of ", feature),
            x = paste0(dim_reduction, "1"),
            y = paste0(dim_reduction, "2"),
            color = "Expression"
          )

        # Add split if requested
        if (!is.null(split_by)) {
          p <- p + ggplot2::facet_wrap(~ .data[[split_by]])
        }

        # Make interactive if requested
        if (interactive) {
          p <- plotly::ggplotly(p)
        }

        # Store plot
        feature_plots[[feature]] <- p
      }

      plot_results$plots$feature_plots <- feature_plots
    }
  }

  # Generate cluster comparison plot
  if ("cluster" %in% plot_types && length(color_by) > 1) {
    # Check for cluster/cell type information
    cluster_cols <- intersect(c("cell_type", "cluster", "seurat_clusters"), colnames(metadata))

    if (length(cluster_cols) > 0) {
      cluster_plots <- list()

      for (cluster_col in cluster_cols) {
        # Calculate contingency table
        if ("modality" %in% colnames(metadata)) {
          contingency_table <- table(
            Cluster = metadata[[cluster_col]],
            Modality = metadata[["modality"]]
          )

          # Convert to data frame for plotting
          table_df <- as.data.frame(contingency_table)
          names(table_df) <- c("Cluster", "Modality", "Count")

          # Create heatmap
          p <- ggplot2::ggplot(table_df, ggplot2::aes(x = Modality, y = Cluster, fill = Count)) +
            ggplot2::geom_tile() +
            ggplot2::scale_fill_viridis_c() +
            ggplot2::theme_minimal() +
            ggplot2::labs(
              title = paste0("Cluster Distribution Across Modalities"),
              x = "Modality",
              y = cluster_col,
              fill = "Cell Count"
            )

          # Make interactive if requested
          if (interactive) {
            p <- plotly::ggplotly(p)
          }

          # Store plot
          cluster_plots[[cluster_col]] <- p
        }
      }

      plot_results$plots$cluster_plots <- cluster_plots
    }
  }

  # Generate density plots
  if ("density" %in% plot_types && "modality" %in% colnames(metadata)) {
    # Generate 2D density plots for each modality
    modalities <- unique(metadata$modality)
    density_plots <- list()

    if (length(modalities) > 1) {
      # Create separate density plots for each modality
      for (mod in modalities) {
        mod_data <- plot_df[plot_df$modality == mod, ]

        # Create density plot
        p <- ggplot2::ggplot(mod_data, ggplot2::aes(x = Dim1, y = Dim2)) +
          ggplot2::stat_density_2d(ggplot2::aes(fill = ..density..), geom = "raster", contour = FALSE) +
          ggplot2::scale_fill_viridis_c() +
          ggplot2::theme_minimal() +
          ggplot2::labs(
            title = paste0("Density Distribution for ", mod),
            x = paste0(dim_reduction, "1"),
            y = paste0(dim_reduction, "2"),
            fill = "Density"
          )

        # Make interactive if requested
        if (interactive) {
          p <- plotly::ggplotly(p)
        }

        # Store plot
        density_plots[[mod]] <- p
      }

      # Create overlay density plot
      p_overlay <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Dim1, y = Dim2, color = modality)) +
        ggplot2::stat_density_2d() +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Density Overlay of Modalities",
          x = paste0(dim_reduction, "1"),
          y = paste0(dim_reduction, "2"),
          color = "Modality"
        )

      if (interactive) {
        p_overlay <- plotly::ggplotly(p_overlay)
      }

      density_plots[["overlay"]] <- p_overlay
      plot_results$plots$density_plots <- density_plots
    }
  }

  # Generate metrics plots if available
  if ("metrics" %in% plot_types && "mixing" %in% names(integration_results)) {
    # Extract mixing scores
    mixing_data <- integration_results$mixing
    plot_results$plots$metrics_plots <- list(mixing = mixing_data$plot)
  }

  # Combine plots if patchwork is available
  if (requireNamespace("patchwork", quietly = TRUE) && !interactive) {
    combined_plots <- list()

    # Combine dimension reduction plots
    if ("dimplots" %in% names(plot_results$plots)) {
      combined_plots$dimplots <- patchwork::wrap_plots(
        plot_results$plots$dimplots,
        ncol = min(ncol, length(plot_results$plots$dimplots))
      ) +
        patchwork::plot_annotation(title = "Dimension Reduction Visualizations")
    }

    # Combine feature plots
    if ("feature_plots" %in% names(plot_results$plots)) {
      combined_plots$feature_plots <- patchwork::wrap_plots(
        plot_results$plots$feature_plots,
        ncol = min(ncol, length(plot_results$plots$feature_plots))
      ) +
        patchwork::plot_annotation(title = "Feature Expression Patterns")
    }

    # Store combined plots
    plot_results$combined_plots <- combined_plots

    # Create master combined plot if requested
    if (length(combined_plots) > 0) {
      master_plot <- patchwork::wrap_plots(
        unlist(combined_plots, recursive = FALSE),
        ncol = 1
      ) +
        patchwork::plot_annotation(
          title = "Multi-Modal Integration Visualization",
          theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold"))
        )

      plot_results$master_plot <- master_plot
    }
  }

  # Return results
  if (return_plots) {
    return(plot_results)
  } else {
    # Print plots and return invisibly
    for (plot_category in names(plot_results$plots)) {
      for (plot_name in names(plot_results$plots[[plot_category]])) {
        print(plot_results$plots[[plot_category]][[plot_name]])
      }
    }

    # Print combined plots if available
    if ("master_plot" %in% names(plot_results)) {
      print(plot_results$master_plot)
    }

    return(invisible(plot_results))
  }
}
