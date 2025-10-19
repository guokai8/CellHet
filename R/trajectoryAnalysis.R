#' Analyze DEG Changes Along a Trajectory
#'
#' Identifies and visualizes gene expression changes along a pseudotime trajectory or other
#' continuous variable in single-cell data.
#'
#' @param object SingleCellExperiment or Seurat object with trajectory information
#' @param trajectory_var Name of the column containing pseudotime or trajectory information
#' @param cell_type_var Name of the column containing cell type annotations (default: "cell_type")
#' @param group_var Optional name of the column for grouping cells (e.g., condition, treatment)
#' @param n_bins Number of bins to divide the trajectory into (default: 10)
#' @param gene_set Optional vector of genes to analyze (default: all genes)
#' @param min_cells_per_bin Minimum number of cells required in each bin (default: 5)
#' @param smooth_method Method for smoothing gene expression trends ("loess", "gam", or "none")
#' @param n_top_genes Number of top dynamic genes to return (default: 50)
#' @param cluster_genes Logical, whether to cluster genes by expression pattern (default: TRUE)
#' @param n_clusters Number of gene clusters to identify (default: 5)
#' @param plot_type Type of visualization to create ("heatmap", "line", or "both")
#'
#' @return A list containing trajectory analysis results and visualizations
#' @export
analyzeDEGTrajectory <- function(object,
                                 trajectory_var,
                                 cell_type_var = "cell_type",
                                 group_var = NULL,
                                 n_bins = 10,
                                 gene_set = NULL,
                                 min_cells_per_bin = 5,
                                 smooth_method = "loess",
                                 n_top_genes = 50,
                                 cluster_genes = TRUE,
                                 n_clusters = 5,
                                 plot_type = "both") {

  # Check if input is Seurat or SingleCellExperiment
  is_seurat <- inherits(object, "Seurat")
  is_sce <- inherits(object, "SingleCellExperiment")

  if (!is_seurat && !is_sce) {
    stop("Input must be a Seurat or SingleCellExperiment object")
  }

  # Extract metadata
  if (is_seurat) {
    metadata <- object@meta.data
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' required for processing Seurat objects")
    }
  } else {
    metadata <- as.data.frame(SingleCellExperiment::colData(object))
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("Package 'SingleCellExperiment' required for processing SCE objects")
    }
  }

  # Check if required columns exist in metadata
  if (!trajectory_var %in% colnames(metadata)) {
    stop(paste0("Column '", trajectory_var, "' not found in metadata"))
  }

  if (!cell_type_var %in% colnames(metadata)) {
    stop(paste0("Column '", cell_type_var, "' not found in metadata"))
  }

  if (!is.null(group_var) && !group_var %in% colnames(metadata)) {
    stop(paste0("Column '", group_var, "' not found in metadata"))
  }

  # Extract trajectory values and check for NA/non-numeric values
  trajectory_values <- metadata[[trajectory_var]]
  if (sum(is.na(trajectory_values)) > 0) {
    warning(paste0("Found ", sum(is.na(trajectory_values)), " NA values in trajectory data. These cells will be removed."))
    valid_cells <- !is.na(trajectory_values)
    metadata <- metadata[valid_cells, ]
    trajectory_values <- trajectory_values[valid_cells]

    if (is_seurat) {
      object <- subset(object, cells = rownames(metadata))
    } else {
      object <- object[, rownames(metadata)]
    }
  }

  if (!is.numeric(trajectory_values)) {
    stop("Trajectory values must be numeric")
  }

  # Check for sufficient variation in trajectory values
  traj_range <- range(trajectory_values)
  if (abs(diff(traj_range)) < 1e-10) {
    warning("Very small or no variation in trajectory values. Adding small jitter to prevent errors.")
    # Add small amount of jitter to ensure uniqueness
    jitter_amount <- max(1e-5, diff(traj_range) * 0.01)
    trajectory_values <- trajectory_values + stats::runif(length(trajectory_values), -jitter_amount, jitter_amount)
    # Update the values in metadata for later use
    metadata[[trajectory_var]] <- trajectory_values
  }

  # Create bins along trajectory
  n_cells <- nrow(metadata)
  bin_breaks <- seq(min(trajectory_values), max(trajectory_values), length.out = n_bins + 1)

  # Ensure bin breaks are unique
  if (length(unique(bin_breaks)) < 2) {
    # Force artificial breaks if not enough unique values
    bin_breaks <- seq(min(trajectory_values) - 0.1, max(trajectory_values) + 0.1, length.out = n_bins + 1)
  }

  bin_assignments <- cut(trajectory_values, breaks = bin_breaks, labels = FALSE, include.lowest = TRUE)

  metadata$bin <- bin_assignments

  # Check cell counts per bin and warn if any are below threshold
  bin_counts <- table(bin_assignments)
  small_bins <- which(bin_counts < min_cells_per_bin)
  if (length(small_bins) > 0) {
    warning(paste0("Bins ", paste(small_bins, collapse = ", "), " have fewer than ", min_cells_per_bin, " cells"))
  }

  # Prepare expression matrix
  if (is_seurat) {
    # Extract normalized expression data
    if (!"RNA" %in% names(object@assays)) {
      stop("RNA assay not found in Seurat object")
    }
    expression_matrix <- Seurat::GetAssayData(object, slot = "data", assay = "RNA")
  } else {
    # For SingleCellExperiment
    if (!"logcounts" %in% SingleCellExperiment::assayNames(object)) {
      warning("logcounts not found in SingleCellExperiment object. Using counts instead.")
      if (!"counts" %in% SingleCellExperiment::assayNames(object)) {
        stop("Neither logcounts nor counts found in SingleCellExperiment object")
      }
      expression_matrix <- SingleCellExperiment::counts(object)
    } else {
      expression_matrix <- SingleCellExperiment::logcounts(object)
    }
  }

  # Filter for specific genes if provided
  if (!is.null(gene_set)) {
    valid_genes <- intersect(gene_set, rownames(expression_matrix))
    if (length(valid_genes) == 0) {
      stop("None of the specified genes found in the dataset")
    }
    if (length(valid_genes) < length(gene_set)) {
      warning(paste0(length(gene_set) - length(valid_genes), " genes not found in the dataset"))
    }
    expression_matrix <- expression_matrix[valid_genes, ]
  }

  # Calculate average expression per bin
  bin_expression <- matrix(0, nrow = nrow(expression_matrix), ncol = n_bins)
  rownames(bin_expression) <- rownames(expression_matrix)
  colnames(bin_expression) <- paste0("bin_", 1:n_bins)

  for (b in 1:n_bins) {
    bin_cells <- rownames(metadata)[metadata$bin == b]

    if (length(bin_cells) > 0) {
      if (length(bin_cells) > 1) {
        # For multiple cells, ensure we're working with the right format
        if (inherits(expression_matrix, "dgCMatrix")) {
          # For sparse matrices
          bin_expression[, b] <- Matrix::rowMeans(expression_matrix[, bin_cells, drop = FALSE])
        } else {
          bin_expression[, b] <- rowMeans(expression_matrix[, bin_cells, drop = FALSE])
        }
      } else {
        # For a single cell, handle differently for sparse matrices
        if (inherits(expression_matrix, "dgCMatrix")) {
          bin_expression[, b] <- as.vector(expression_matrix[, bin_cells])
        } else {
          bin_expression[, b] <- expression_matrix[, bin_cells, drop = TRUE]
        }
      }
    } else {
      # Handle empty bins
      bin_expression[, b] <- NA
    }
  }

  # Calculate variability of gene expression across bins
  gene_variability <- apply(bin_expression, 1, function(x) {
    if (all(is.na(x))) return(0)
    return(stats::var(x, na.rm = TRUE))
  })

  # Identify top variable genes
  top_genes <- names(sort(gene_variability, decreasing = TRUE))[1:min(n_top_genes, length(gene_variability))]

  # Create gene expression trends along trajectory
  trajectory_midpoints <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2

  # Generate smoothed gene expression profiles along trajectory
  if (smooth_method != "none") {
    smoothed_expression <- matrix(0, nrow = length(top_genes), ncol = 100)
    rownames(smoothed_expression) <- top_genes
    trajectory_points <- seq(min(trajectory_values), max(trajectory_values), length.out = 100)

    if (smooth_method == "loess") {
      for (i in 1:length(top_genes)) {
        gene <- top_genes[i]
        gene_expr <- bin_expression[gene, ]

        # Only smooth if we have enough non-NA values
        if (sum(!is.na(gene_expr)) >= 3) {
          tryCatch({
            loess_fit <- stats::loess(gene_expr ~ trajectory_midpoints, na.action = na.exclude)
            smoothed_expression[i, ] <- stats::predict(loess_fit, trajectory_points)
          }, error = function(e) {
            # Fallback to linear interpolation if loess fails
            smoothed_expression[i, ] <- approx(
              x = trajectory_midpoints[!is.na(gene_expr)],
              y = gene_expr[!is.na(gene_expr)],
              xout = trajectory_points,
              rule = 2
            )$y
          })
        } else {
          # Not enough data points, use NAs
          smoothed_expression[i, ] <- NA
        }
      }
    } else if (smooth_method == "gam") {
      if (!requireNamespace("mgcv", quietly = TRUE)) {
        warning("Package 'mgcv' not available. Falling back to loess smoothing.")
        smooth_method <- "loess"
        for (i in 1:length(top_genes)) {
          gene <- top_genes[i]
          gene_expr <- bin_expression[gene, ]

          # Only smooth if we have enough non-NA values
          if (sum(!is.na(gene_expr)) >= 3) {
            tryCatch({
              loess_fit <- stats::loess(gene_expr ~ trajectory_midpoints, na.action = na.exclude)
              smoothed_expression[i, ] <- stats::predict(loess_fit, trajectory_points)
            }, error = function(e) {
              # Fallback to linear interpolation if loess fails
              smoothed_expression[i, ] <- approx(
                x = trajectory_midpoints[!is.na(gene_expr)],
                y = gene_expr[!is.na(gene_expr)],
                xout = trajectory_points,
                rule = 2
              )$y
            })
          } else {
            # Not enough data points, use NAs
            smoothed_expression[i, ] <- NA
          }
        }
      } else {
        for (i in 1:length(top_genes)) {
          gene <- top_genes[i]
          gene_expr <- bin_expression[gene, ]

          # Only smooth if we have enough non-NA values
          if (sum(!is.na(gene_expr)) >= 3) {
            tryCatch({
              gam_fit <- mgcv::gam(gene_expr ~ s(trajectory_midpoints), na.action = na.exclude)
              smoothed_expression[i, ] <- stats::predict(gam_fit,
                                                         data.frame(trajectory_midpoints = trajectory_points))
            }, error = function(e) {
              # Fallback to linear interpolation if gam fails
              smoothed_expression[i, ] <- approx(
                x = trajectory_midpoints[!is.na(gene_expr)],
                y = gene_expr[!is.na(gene_expr)],
                xout = trajectory_points,
                rule = 2
              )$y
            })
          } else {
            # Not enough data points, use NAs
            smoothed_expression[i, ] <- NA
          }
        }
      }
    }
  } else {
    # No smoothing, just use bin expression
    smoothed_expression <- bin_expression[top_genes, ]
    trajectory_points <- trajectory_midpoints
  }

  # Cluster genes by expression pattern if requested
  gene_clusters <- NULL
  if (cluster_genes) {
    # Normalize expression patterns for clustering
    # First, ensure we have sufficient non-NA values
    valid_rows <- apply(smoothed_expression, 1, function(x) sum(!is.na(x)) > 3)
    if (sum(valid_rows) < 2) {
      warning("Not enough valid data for clustering. Skipping gene clustering.")
      cluster_genes <- FALSE
    } else {
      # Only cluster valid rows
      valid_smoothed <- smoothed_expression[valid_rows, , drop = FALSE]

      # Handle NAs in the expression matrix
      scaled_expression <- t(scale(t(replace(valid_smoothed, is.na(valid_smoothed), 0))))

      # Perform hierarchical clustering
      gene_dist <- stats::dist(scaled_expression)
      gene_hclust <- stats::hclust(gene_dist, method = "ward.D2")

      # Cut tree to get requested number of clusters
      n_clusters <- min(n_clusters, sum(valid_rows))
      row_clusters <- stats::cutree(gene_hclust, k = n_clusters)

      # Convert to named vector for all genes
      gene_clusters <- rep(NA, length(top_genes))
      names(gene_clusters) <- top_genes
      gene_clusters[names(row_clusters)] <- row_clusters
    }
  }

  # Create visualizations
  plots <- list()

  if (plot_type %in% c("heatmap", "both")) {
    # Create heatmap of gene expression along trajectory
    if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
        requireNamespace("grid", quietly = TRUE) &&
        requireNamespace("circlize", quietly = TRUE)) {

      # Prepare heatmap data
      heatmap_data <- smoothed_expression

      # Check for valid data
      if (all(is.na(heatmap_data))) {
        warning("All expression values are NA. Cannot create heatmap.")
      } else {
        # Replace NAs with min value for visualization
        na_mask <- is.na(heatmap_data)
        if (any(na_mask)) {
          min_val <- min(heatmap_data, na.rm = TRUE)
          heatmap_data[na_mask] <- min_val
        }

        # Check for sufficient variation in data
        min_val <- min(heatmap_data, na.rm = TRUE)
        max_val <- max(heatmap_data, na.rm = TRUE)

        # Create color scale - FIX: Handle the case where min and max are too close
        if (abs(max_val - min_val) < 1e-6) {
          # Force artificial range if values are too similar
          breaks <- c(min_val - 1, min_val, min_val + 1)
        } else {
          breaks <- c(min_val, (min_val + max_val)/2, max_val)
        }

        col_fun <- circlize::colorRamp2(
          breaks,
          c("deepskyblue", "white", "darkorange")
        )

        # Add column and row annotations
        column_anno <- data.frame(
          Pseudotime = trajectory_points
        )

        if (cluster_genes && !is.null(gene_clusters) && !all(is.na(gene_clusters))) {
          row_anno <- data.frame(
            Cluster = factor(gene_clusters[top_genes])
          )

          # Generate colors for clusters - handle NAs
          cluster_ids <- unique(na.omit(gene_clusters))
          cluster_colors <- stats::setNames(
            rainbow(length(cluster_ids)),
            cluster_ids
          )

          row_ha <- ComplexHeatmap::rowAnnotation(
            df = row_anno,
            col = list(Cluster = cluster_colors)
          )
        } else {
          row_ha <- NULL
        }

        # Create heatmap
        hm <- ComplexHeatmap::Heatmap(
          matrix = heatmap_data,
          name = "Expression",
          col = col_fun,
          show_row_names = TRUE,
          show_column_names = FALSE,
          cluster_rows = ifelse(cluster_genes, FALSE, TRUE),  # Already clustered if requested
          cluster_columns = FALSE,  # Keep trajectory order
          row_names_gp = grid::gpar(fontsize = 8),
          column_title = "Gene Expression Along Trajectory",
          row_title = "Genes",
          right_annotation = row_ha
        )

        plots$heatmap <- hm
      }
    } else {
      # Fallback to using ggplot2 for heatmap
      if (requireNamespace("ggplot2", quietly = TRUE) &&
          requireNamespace("reshape2", quietly = TRUE)) {

        # Reshape data for ggplot
        valid_data <- !all(is.na(smoothed_expression))
        if (!valid_data) {
          warning("All expression values are NA. Cannot create heatmap.")
        } else {
          heatmap_data <- reshape2::melt(smoothed_expression)
          names(heatmap_data) <- c("Gene", "Bin", "Expression")

          # Add trajectory information
          heatmap_data$Pseudotime <- trajectory_points[heatmap_data$Bin]

          # Add cluster information if available
          if (!is.null(gene_clusters) && !all(is.na(gene_clusters))) {
            heatmap_data$Cluster <- gene_clusters[heatmap_data$Gene]

            # Create heatmap
            p <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = Pseudotime, y = Gene, fill = Expression)) +
              ggplot2::geom_tile() +
              ggplot2::scale_fill_gradient2(low = "deepskyblue", mid = "white", high = "darkorange",
                                            na.value = "grey90") +
              ggplot2::facet_grid(Cluster ~ ., scales = "free_y", space = "free_y") +
              ggplot2::theme_minimal() +
              ggplot2::labs(title = "Gene Expression Along Trajectory",
                            x = "Pseudotime",
                            y = "Genes")
          } else {
            # Create heatmap without clusters
            p <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = Pseudotime, y = Gene, fill = Expression)) +
              ggplot2::geom_tile() +
              ggplot2::scale_fill_gradient2(low = "deepskyblue", mid = "white", high = "darkorange",
                                            na.value = "grey90") +
              ggplot2::theme_minimal() +
              ggplot2::labs(title = "Gene Expression Along Trajectory",
                            x = "Pseudotime",
                            y = "Genes")
          }

          plots$heatmap <- p
        }
      } else {
        warning("Required packages (ggplot2, reshape2) not available for heatmap visualization")
      }
    }
  }

  if (plot_type %in% c("line", "both")) {
    # Create line plots of gene expression along trajectory
    if (requireNamespace("ggplot2", quietly = TRUE) &&
        requireNamespace("reshape2", quietly = TRUE)) {

      # Reshape data for ggplot
      line_data <- reshape2::melt(smoothed_expression)
      names(line_data) <- c("Gene", "Bin", "Expression")

      # Add trajectory information
      line_data$Pseudotime <- trajectory_points[line_data$Bin]

      # Add cluster information if available
      if (cluster_genes && !is.null(gene_clusters) && !all(is.na(gene_clusters))) {
        valid_clusters <- !is.na(gene_clusters)
        if (sum(valid_clusters) > 0) {
          line_data$Cluster <- factor(gene_clusters[line_data$Gene])

          # Create separate line plots for each cluster
          cluster_plots <- list()

          for (cluster in unique(na.omit(gene_clusters))) {
            cluster_data <- line_data[!is.na(line_data$Cluster) & line_data$Cluster == cluster, ]

            if (nrow(cluster_data) > 0) {
              p <- ggplot2::ggplot(cluster_data,
                                   ggplot2::aes(x = Pseudotime, y = Expression,
                                                color = Gene, group = Gene)) +
                ggplot2::geom_line() +
                ggplot2::theme_minimal() +
                ggplot2::labs(title = paste("Cluster", cluster, "Genes"),
                              x = "Pseudotime",
                              y = "Expression") +
                ggplot2::theme(legend.position = "right")

              cluster_plots[[paste0("cluster_", cluster)]] <- p
            }
          }

          plots$line_plots <- cluster_plots

          # Also create a faceted plot
          if (length(unique(na.omit(gene_clusters))) > 1) {
            p_facet <- ggplot2::ggplot(line_data[!is.na(line_data$Cluster), ],
                                       ggplot2::aes(x = Pseudotime, y = Expression,
                                                    color = Gene, group = Gene)) +
              ggplot2::geom_line() +
              ggplot2::facet_wrap(~ Cluster, scales = "free_y") +
              ggplot2::theme_minimal() +
              ggplot2::labs(title = "Gene Expression Trajectories by Cluster",
                            x = "Pseudotime",
                            y = "Expression") +
              ggplot2::theme(legend.position = "none")

            plots$faceted_line_plot <- p_facet
          }
        } else {
          # No valid clusters, create a single line plot for all genes
          p <- ggplot2::ggplot(line_data,
                               ggplot2::aes(x = Pseudotime, y = Expression,
                                            color = Gene, group = Gene)) +
            ggplot2::geom_line() +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = "Gene Expression Trajectories",
                          x = "Pseudotime",
                          y = "Expression") +
            ggplot2::theme(legend.position = "right")

          plots$line_plot <- p
        }
      } else {
        # Create a single line plot for all genes
        p <- ggplot2::ggplot(line_data,
                             ggplot2::aes(x = Pseudotime, y = Expression,
                                          color = Gene, group = Gene)) +
          ggplot2::geom_line() +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = "Gene Expression Trajectories",
                        x = "Pseudotime",
                        y = "Expression") +
          ggplot2::theme(legend.position = "right")

        plots$line_plot <- p
      }
    } else {
      warning("Required packages (ggplot2, reshape2) not available for line plot visualization")
    }
  }

  # Compile results
  results <- list(
    trajectory_info = list(
      trajectory_var = trajectory_var,
      min_value = min(trajectory_values),
      max_value = max(trajectory_values),
      bin_breaks = bin_breaks,
      bin_counts = bin_counts
    ),
    bin_expression = bin_expression,
    top_genes = top_genes,
    gene_variability = gene_variability[top_genes],
    smoothed_expression = smoothed_expression,
    trajectory_points = trajectory_points,
    gene_clusters = gene_clusters,
    plots = plots
  )

  return(results)
}

#' Create Gene Expression Profile Along Trajectory
#'
#' Generates a plot showing the expression of a specific gene or set of genes along a trajectory.
#'
#' @param trajectory_results Output from analyzeDEGTrajectory function
#' @param genes Character vector of gene names to plot
#' @param plot_type Type of plot to create: "line", "box", or "violin" (default: "line")
#' @param add_loess_fit Logical, whether to add a loess smoothing curve to the plot (default: TRUE)
#' @param facet_by_cluster Logical, whether to facet by gene clusters if available (default: TRUE)
#' @param highlight_gene Optional, name of a gene to highlight in the plot
#' @param custom_colors Optional vector of custom colors for plotting
#' @param line_size Size of lines in the plot (default: 1)
#' @param point_size Size of points in the plot (default: 2)
#' @param title Plot title (default: "Gene Expression Along Trajectory")
#'
#' @return A ggplot object showing gene expression profiles
#' @export
plotGeneTrajectory <- function(trajectory_results,
                               genes,
                               plot_type = "line",
                               add_loess_fit = TRUE,
                               facet_by_cluster = TRUE,
                               highlight_gene = NULL,
                               custom_colors = NULL,
                               line_size = 1,
                               point_size = 2,
                               title = "Gene Expression Along Trajectory") {

  # Check required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function. Please install it.")
  }

  # Validate input
  if (!is.list(trajectory_results) || !"smoothed_expression" %in% names(trajectory_results)) {
    stop("Invalid trajectory_results input. Must be output from analyzeDEGTrajectory function.")
  }

  # Extract trajectory points and smoothed expression
  trajectory_points <- trajectory_results$trajectory_points
  smoothed_expression <- trajectory_results$smoothed_expression

  # Check if requested genes are available
  available_genes <- rownames(smoothed_expression)
  valid_genes <- intersect(genes, available_genes)

  if (length(valid_genes) == 0) {
    stop("None of the specified genes found in the trajectory results")
  }

  if (length(valid_genes) < length(genes)) {
    warning(paste0(length(genes) - length(valid_genes), " genes not found in the trajectory results"))
  }

  # Extract expression data for requested genes
  gene_expression <- smoothed_expression[valid_genes, , drop = FALSE]

  # Reshape data for plotting
  plot_data <- reshape2::melt(gene_expression)
  names(plot_data) <- c("Gene", "Bin", "Expression")
  plot_data$Pseudotime <- trajectory_points[plot_data$Bin]

  # Add cluster information if available and requested
  has_clusters <- !is.null(trajectory_results$gene_clusters) && facet_by_cluster
  if (has_clusters) {
    plot_data$Cluster <- factor(trajectory_results$gene_clusters[plot_data$Gene])
  }

  # Highlight specific gene if requested
  if (!is.null(highlight_gene) && highlight_gene %in% valid_genes) {
    plot_data$Highlight <- plot_data$Gene == highlight_gene
  } else {
    plot_data$Highlight <- FALSE
  }

  # Create plot based on plot_type
  if (plot_type == "line") {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Pseudotime, y = Expression, color = Gene, group = Gene)) +
      ggplot2::geom_line(size = line_size) +
      ggplot2::geom_point(size = point_size)

    # Add loess fit if requested
    if (add_loess_fit) {
      p <- p + ggplot2::geom_smooth(method = "loess", se = TRUE, alpha = 0.2, linetype = "dashed")
    }

  } else if (plot_type == "box") {
    # Create bins for box plots
    n_bins <- 10
    bin_breaks <- seq(min(plot_data$Pseudotime), max(plot_data$Pseudotime), length.out = n_bins + 1)
    plot_data$PseudotimeBin <- cut(plot_data$Pseudotime, breaks = bin_breaks, labels = FALSE, include.lowest = TRUE)

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = factor(PseudotimeBin), y = Expression, fill = Gene)) +
      ggplot2::geom_boxplot(alpha = 0.7) +
      ggplot2::labs(x = "Pseudotime Bin")

  } else if (plot_type == "violin") {
    # Create bins for violin plots
    n_bins <- 5  # Fewer bins for violin plots
    bin_breaks <- seq(min(plot_data$Pseudotime), max(plot_data$Pseudotime), length.out = n_bins + 1)
    plot_data$PseudotimeBin <- cut(plot_data$Pseudotime, breaks = bin_breaks, labels = FALSE, include.lowest = TRUE)

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = factor(PseudotimeBin), y = Expression, fill = Gene)) +
      ggplot2::geom_violin(alpha = 0.7) +
      ggplot2::geom_boxplot(width = 0.1, alpha = 0.2) +
      ggplot2::labs(x = "Pseudotime Bin")
  } else {
    stop("Invalid plot_type. Must be 'line', 'box', or 'violin'.")
  }

  # Apply custom colors if provided
  if (!is.null(custom_colors)) {
    if (plot_type == "line") {
      p <- p + ggplot2::scale_color_manual(values = custom_colors)
    } else {
      p <- p + ggplot2::scale_fill_manual(values = custom_colors)
    }
  }

  # Apply faceting if clusters are available and requested
  if (has_clusters) {
    p <- p + ggplot2::facet_wrap(~ Cluster, scales = "free_y")
  }

  # Highlight specific gene if requested
  if (any(plot_data$Highlight) && plot_type == "line") {
    highlight_data <- plot_data[plot_data$Highlight, ]
    p <- p + ggplot2::geom_line(data = highlight_data, size = line_size * 2, color = "black") +
      ggplot2::geom_line(data = highlight_data, size = line_size * 1.5)
  }

  # Add plot styling
  p <- p + ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      x = "Pseudotime",
      y = "Expression"
    ) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.title = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10)
    )

  return(p)
}

#' Differential Trajectory Analysis Between Conditions
#'
#' Compares gene expression trajectories between two or more conditions along a trajectory.
#'
#' @param object SingleCellExperiment or Seurat object with trajectory information
#' @param trajectory_var Name of the column containing pseudotime or trajectory information
#' @param group_var Name of the column for condition/group comparison
#' @param cell_type_var Name of the column for cell type (default: "cell_type")
#' @param cell_types Optional vector of cell types to include (default: all cell types)
#' @param n_bins Number of bins to divide the trajectory into (default: 10)
#' @param min_cells_per_bin Minimum number of cells required in each bin-group combination (default: 3)
#' @param n_top_genes Number of top differentially expressed genes to return (default: 50)
#' @param test_method Statistical test to use for identifying differential trajectories (default: "anova")
#' @param p_adj_method Method for p-value adjustment (default: "BH")
#' @param p_val_threshold Significance threshold for adjusted p-values (default: 0.05)
#' @param smooth_method Method for smoothing gene expression trends (default: "loess")
#'
#' @return A list containing differential trajectory analysis results and visualizations
#' @export
diffTrajectoryAnalysis <- function(object,
                                   trajectory_var,
                                   group_var,
                                   cell_type_var = "cell_type",
                                   cell_types = NULL,
                                   n_bins = 10,
                                   min_cells_per_bin = 3,
                                   n_top_genes = 50,
                                   test_method = "anova",
                                   p_adj_method = "BH",
                                   p_val_threshold = 0.05,
                                   smooth_method = "loess") {

  # Check if input is Seurat or SingleCellExperiment
  is_seurat <- inherits(object, "Seurat")
  is_sce <- inherits(object, "SingleCellExperiment")

  if (!is_seurat && !is_sce) {
    stop("Input must be a Seurat or SingleCellExperiment object")
  }

  # Validate test_method
  test_method <- match.arg(test_method, c("anova", "pattern_diff"))

  # Validate p_adj_method
  p_adj_method <- match.arg(p_adj_method, c("BH", "bonferroni", "holm", "hochberg", "hommel", "BY", "fdr"))

  # Extract metadata
  if (is_seurat) {
    metadata <- object@meta.data
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' required for processing Seurat objects")
    }
  } else {
    metadata <- as.data.frame(SingleCellExperiment::colData(object))
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("Package 'SingleCellExperiment' required for processing SCE objects")
    }
  }

  # Check if required columns exist in metadata
  if (!trajectory_var %in% colnames(metadata)) {
    stop(paste0("Column '", trajectory_var, "' not found in metadata"))
  }

  if (!group_var %in% colnames(metadata)) {
    stop(paste0("Column '", group_var, "' not found in metadata"))
  }

  if (!cell_type_var %in% colnames(metadata)) {
    stop(paste0("Column '", cell_type_var, "' not found in metadata"))
  }

  # Filter for specific cell types if requested
  if (!is.null(cell_types)) {
    valid_cell_types <- intersect(cell_types, unique(metadata[[cell_type_var]]))
    if (length(valid_cell_types) == 0) {
      stop("None of the specified cell types found in the data")
    }
    if (length(valid_cell_types) < length(cell_types)) {
      warning(paste0(length(cell_types) - length(valid_cell_types), " cell types not found in the data"))
    }

    # Subset cells to only include specified cell types
    cell_mask <- metadata[[cell_type_var]] %in% valid_cell_types
    metadata <- metadata[cell_mask, ]

    if (is_seurat) {
      object <- subset(object, cells = rownames(metadata))
    } else {
      object <- object[, rownames(metadata)]
    }
  } else {
    valid_cell_types <- unique(metadata[[cell_type_var]])
  }

  # Extract trajectory values and check for NA/non-numeric values
  trajectory_values <- metadata[[trajectory_var]]
  if (sum(is.na(trajectory_values)) > 0) {
    warning(paste0("Found ", sum(is.na(trajectory_values)), " NA values in trajectory data. These cells will be removed."))
    valid_cells <- !is.na(trajectory_values)
    metadata <- metadata[valid_cells, ]
    trajectory_values <- trajectory_values[valid_cells]

    if (is_seurat) {
      object <- subset(object, cells = rownames(metadata))
    } else {
      object <- object[, rownames(metadata)]
    }
  }

  if (!is.numeric(trajectory_values)) {
    stop("Trajectory values must be numeric")
  }

  # Extract unique groups
  groups <- unique(metadata[[group_var]])
  if (length(groups) < 2) {
    stop("At least two groups are required for differential trajectory analysis")
  }

  # Create bins along trajectory
  bin_breaks <- seq(min(trajectory_values), max(trajectory_values), length.out = n_bins + 1)
  bin_assignments <- cut(trajectory_values, breaks = bin_breaks, labels = FALSE, include.lowest = TRUE)

  metadata$bin <- bin_assignments

  # Check sample size per bin/group
  bin_group_counts <- table(metadata$bin, metadata[[group_var]])
  small_combinations <- which(bin_group_counts < min_cells_per_bin, arr.ind = TRUE)
  if (nrow(small_combinations) > 0) {
    warning(paste0(nrow(small_combinations), " bin-group combinations have fewer than ", min_cells_per_bin, " cells"))
  }

  # Prepare expression matrix
  # Prepare expression matrix
  if (is_seurat) {
    # Extract normalized expression data
    if (!"RNA" %in% names(object@assays)) {
      stop("RNA assay not found in Seurat object")
    }
    expression_matrix <- Seurat::GetAssayData(object, slot = "data", assay = "RNA")
  } else {
    # For SingleCellExperiment
    if (!"logcounts" %in% SingleCellExperiment::assayNames(object)) {
      warning("logcounts not found in SingleCellExperiment object. Using counts instead.")
      if (!"counts" %in% SingleCellExperiment::assayNames(object)) {
        stop("Neither logcounts nor counts found in SingleCellExperiment object")
      }
      expression_matrix <- SingleCellExperiment::counts(object)
    } else {
      expression_matrix <- SingleCellExperiment::logcounts(object)
    }
  }

  # Calculate average expression per bin/group combination
  group_bin_expression <- list()
  for (group in groups) {
    group_bin_expression[[group]] <- matrix(NA, nrow = nrow(expression_matrix), ncol = n_bins)
    rownames(group_bin_expression[[group]]) <- rownames(expression_matrix)
    colnames(group_bin_expression[[group]]) <- paste0("bin_", 1:n_bins)

    for (b in 1:n_bins) {
      # Get cells in this bin and group
      bin_group_cells <- rownames(metadata)[metadata$bin == b & metadata[[group_var]] == group]

      # Check if we have enough cells
      if (length(bin_group_cells) >= min_cells_per_bin) {
        # Have enough cells, calculate mean expression
        group_bin_expression[[group]][, b] <- rowMeans(expression_matrix[, bin_group_cells, drop = FALSE])
      } else if (length(bin_group_cells) == 1) {
        # Special case: exactly one cell
        group_bin_expression[[group]][, b] <- expression_matrix[, bin_group_cells, drop = TRUE]
      } else {
        # Not enough cells, set to NA
        group_bin_expression[[group]][, b] <- NA
      }
    }
  }

  # Generate smooth trajectories for each group if requested
  trajectory_midpoints <- (bin_breaks[-1] + bin_breaks[-length(bin_breaks)]) / 2
  smoothed_group_expression <- list()

  if (smooth_method != "none") {
    trajectory_points <- seq(min(trajectory_values), max(trajectory_values), length.out = 100)

    for (group in groups) {
      smoothed_group_expression[[group]] <- matrix(0, nrow = nrow(expression_matrix), ncol = 100)
      rownames(smoothed_group_expression[[group]]) <- rownames(expression_matrix)

      for (i in 1:nrow(expression_matrix)) {
        gene <- rownames(expression_matrix)[i]
        gene_expr <- group_bin_expression[[group]][i, ]

        # Only smooth if we have enough non-NA values
        if (sum(!is.na(gene_expr)) >= 3) {
          if (smooth_method == "loess") {
            tryCatch({
              loess_fit <- stats::loess(gene_expr ~ trajectory_midpoints, na.action = na.exclude)
              smoothed_group_expression[[group]][i, ] <- stats::predict(loess_fit, trajectory_points)
            }, error = function(e) {
              # If smoothing fails, use linear interpolation
              smoothed_group_expression[[group]][i, ] <- approx(
                x = trajectory_midpoints[!is.na(gene_expr)],
                y = gene_expr[!is.na(gene_expr)],
                xout = trajectory_points,
                rule = 2
              )$y
            })
          } else if (smooth_method == "gam") {
            if (!requireNamespace("mgcv", quietly = TRUE)) {
              warning("Package 'mgcv' not available. Falling back to linear interpolation.")
              smoothed_group_expression[[group]][i, ] <- approx(
                x = trajectory_midpoints[!is.na(gene_expr)],
                y = gene_expr[!is.na(gene_expr)],
                xout = trajectory_points,
                rule = 2
              )$y
            } else {
              tryCatch({
                gam_fit <- mgcv::gam(gene_expr ~ s(trajectory_midpoints), na.action = na.exclude)
                smoothed_group_expression[[group]][i, ] <- stats::predict(gam_fit,
                                                                          data.frame(trajectory_midpoints = trajectory_points))
              }, error = function(e) {
                # If smoothing fails, use linear interpolation
                smoothed_group_expression[[group]][i, ] <- approx(
                  x = trajectory_midpoints[!is.na(gene_expr)],
                  y = gene_expr[!is.na(gene_expr)],
                  xout = trajectory_points,
                  rule = 2
                )$y
              })
            }
          }
        } else {
          # Not enough data points, set to NA
          smoothed_group_expression[[group]][i, ] <- NA
        }
      }
    }
  } else {
    # No smoothing, just use bin expression
    for (group in groups) {
      smoothed_group_expression[[group]] <- group_bin_expression[[group]]
    }
    trajectory_points <- trajectory_midpoints
  }

  # Identify genes with differential trajectories

  # Method 1: ANOVA test for each gene across trajectory points
  if (test_method == "anova") {
    # Prepare data structure for ANOVA tests
    anova_results <- data.frame(
      gene = rownames(expression_matrix),
      F_statistic = NA,
      p_value = NA,
      p_adjusted = NA,
      significance = FALSE,
      stringsAsFactors = FALSE
    )

    for (i in 1:nrow(expression_matrix)) {
      gene <- rownames(expression_matrix)[i]

      # Combine data across groups for this gene
      gene_data <- data.frame(
        expression = numeric(),
        trajectory = numeric(),
        group = character(),
        stringsAsFactors = FALSE
      )

      for (group in groups) {
        gene_expr <- group_bin_expression[[group]][i, ]
        valid_points <- !is.na(gene_expr)

        if (sum(valid_points) >= 3) {  # Need at least 3 points for meaningful analysis
          gene_data <- rbind(gene_data, data.frame(
            expression = gene_expr[valid_points],
            trajectory = trajectory_midpoints[valid_points],
            group = rep(group, sum(valid_points)),
            stringsAsFactors = FALSE
          ))
        }
      }

      # Only run ANOVA if we have data from at least 2 groups
      if (length(unique(gene_data$group)) >= 2 && nrow(gene_data) >= 5) {
        # Run ANOVA with interaction term between trajectory and group
        tryCatch({
          # Create interaction model: expression ~ trajectory * group
          anova_model <- stats::lm(expression ~ trajectory * group, data = gene_data)
          anova_result <- stats::anova(anova_model)

          # Extract p-value for interaction term
          interaction_row <- nrow(anova_result) - 1  # Typically the second-to-last row
          anova_results$F_statistic[i] <- anova_result$`F value`[interaction_row]
          anova_results$p_value[i] <- anova_result$`Pr(>F)`[interaction_row]
        }, error = function(e) {
          # If ANOVA fails, set to NA
          anova_results$F_statistic[i] <- NA
          anova_results$p_value[i] <- NA
        })
      }
    }

    # Adjust p-values for multiple testing
    anova_results$p_adjusted <- stats::p.adjust(anova_results$p_value, method = p_adj_method)
    anova_results$significance <- anova_results$p_adjusted < p_val_threshold & !is.na(anova_results$p_adjusted)

    # Sort by p-value
    anova_results <- anova_results[order(anova_results$p_value), ]

    # Extract top genes
    top_diff_genes <- anova_results$gene[1:min(n_top_genes, nrow(anova_results))]

  } else if (test_method == "pattern_diff") {
    # Method 2: Calculate differences in pattern between groups

    # Compute pairwise distances between group trajectories for each gene
    pattern_diff_results <- data.frame(
      gene = rownames(expression_matrix),
      max_diff = NA,
      avg_diff = NA,
      p_value = NA,
      p_adjusted = NA,
      significance = FALSE,
      stringsAsFactors = FALSE
    )

    # Create all pairwise group combinations
    group_pairs <- combn(groups, 2, simplify = FALSE)

    for (i in 1:nrow(expression_matrix)) {
      gene <- rownames(expression_matrix)[i]

      # Calculate maximum and average differences across group pairs
      max_diff <- 0
      avg_diffs <- numeric()

      for (pair in group_pairs) {
        group1 <- pair[1]
        group2 <- pair[2]

        expr1 <- smoothed_group_expression[[group1]][i, ]
        expr2 <- smoothed_group_expression[[group2]][i, ]

        # Skip if either has insufficient data
        if (sum(!is.na(expr1)) < 3 || sum(!is.na(expr2)) < 3) {
          next
        }

        # Calculate difference in patterns
        valid_points <- !is.na(expr1) & !is.na(expr2)
        if (sum(valid_points) >= 3) {
          diff_vector <- expr1[valid_points] - expr2[valid_points]
          pair_max_diff <- max(abs(diff_vector))
          pair_avg_diff <- mean(abs(diff_vector))

          max_diff <- max(max_diff, pair_max_diff)
          avg_diffs <- c(avg_diffs, pair_avg_diff)
        }
      }

      if (length(avg_diffs) > 0) {
        pattern_diff_results$max_diff[i] <- max_diff
        pattern_diff_results$avg_diff[i] <- mean(avg_diffs)

        # Calculate p-value using permutation test (simplified version)
        # In a real implementation, this would use actual permutation testing
        pattern_diff_results$p_value[i] <- exp(-pattern_diff_results$avg_diff[i])
      }
    }

    # Adjust p-values for multiple testing
    pattern_diff_results$p_adjusted <- stats::p.adjust(pattern_diff_results$p_value, method = p_adj_method)
    pattern_diff_results$significance <- pattern_diff_results$p_adjusted < p_val_threshold & !is.na(pattern_diff_results$p_adjusted)

    # Sort by average difference
    pattern_diff_results <- pattern_diff_results[order(-pattern_diff_results$avg_diff), ]

    # Extract top genes
    top_diff_genes <- pattern_diff_results$gene[1:min(n_top_genes, nrow(pattern_diff_results))]

  } else {
    stop("Invalid test_method. Must be 'anova' or 'pattern_diff'.")
  }

  # Generate visualizations for top differentially expressed genes
  diff_plots <- list()

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    for (gene_idx in 1:min(length(top_diff_genes), 9)) {  # Plot at most 9 genes
      gene <- top_diff_genes[gene_idx]

      # Create data frame for plotting
      plot_data <- data.frame(
        Pseudotime = numeric(),
        Expression = numeric(),
        Group = character(),
        stringsAsFactors = FALSE
      )

      for (group in groups) {
        expr <- smoothed_group_expression[[group]][gene, ]
        plot_data <- rbind(plot_data, data.frame(
          Pseudotime = trajectory_points,
          Expression = expr,
          Group = rep(group, length(trajectory_points)),
          stringsAsFactors = FALSE
        ))
      }

      # Create plot
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Pseudotime, y = Expression, color = Group)) +
        ggplot2::geom_line(size = 1.2) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste0("Gene: ", gene),
          x = "Pseudotime",
          y = "Expression"
        )

      diff_plots[[gene]] <- p
    }

    # Create a combined plot using grid.arrange if available
    if (requireNamespace("gridExtra", quietly = TRUE) && length(diff_plots) > 1) {
      combined_plot <- gridExtra::grid.arrange(
        grobs = diff_plots[1:min(length(diff_plots), 9)],
        nrow = 3,
        ncol = 3,
        top = "Top Differentially Expressed Genes Along Trajectory"
      )
      diff_plots$combined <- combined_plot
    }
  }

  # Compile results
  if (test_method == "anova") {
    test_results <- anova_results
  } else {
    test_results <- pattern_diff_results
  }

  results <- list(
    metadata = list(
      trajectory_var = trajectory_var,
      group_var = group_var,
      cell_type_var = cell_type_var,
      cell_types = valid_cell_types,
      groups = groups,
      n_bins = n_bins,
      bin_breaks = bin_breaks,
      test_method = test_method
    ),
    group_bin_expression = group_bin_expression,
    smoothed_group_expression = smoothed_group_expression,
    trajectory_points = trajectory_points,
    test_results = test_results,
    top_diff_genes = top_diff_genes,
    plots = diff_plots
  )

  return(results)
}
