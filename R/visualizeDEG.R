#' Create a DEG Heatmap for Cell Types and Comparisons
#'
#' @param deg_results Output from findDifferentialGenes or compareDEGs function
#' @param comparison_col Column name for the comparison categories
#' @param cell_type_col Column name for the cell types
#' @param border_color Color for tile borders (default = "white")
#' @param border_size Size of tile borders (default = 0.5)
#' @param show_counts Logical. Whether to show count numbers on tiles (default = TRUE)
#' @param count_size Size of count numbers (default = 3)
#' @param count_color Color of count numbers (default = "black")
#' @param low_color Color for low values (default = "deepskyblue")
#' @param mid_color Color for mid values (default = "white")
#' @param high_color Color for high values (default = "darkorange")
#' @param cluster_rows Logical. Whether to cluster rows (default = TRUE)
#' @param cluster_cols Logical. Whether to cluster columns (default = FALSE)
#' @param cluster_method Clustering method to use (default = "complete")
#' @param distance_method Distance method for clustering (default = "euclidean")
#' @param group_regulation Logical. If TRUE, groups up/down on x-axis (default = FALSE)
#' @param show_reg_indicators Logical. If TRUE, adds color indicators for up/down regulation (default = TRUE)
#' @param up_indicator_color Color for up-regulation indicator (default = "darkorange")
#' @param down_indicator_color Color for down-regulation indicator (default = "deepskyblue")
#' @param indicator_height Height of the color indicator (default = 0.2)
#' @param top_labels Logical. If TRUE, shows Up/Down labels at top of heatmap (default = TRUE)
#' @param top_label_size Size of the top labels (default = 3.5)
#' @param label_distance Distance of labels from the top of the heatmap (default = 1.5)
#' @param title Plot title (default = "DEGs Across Cell Type and Comparisons")
#' @param title_size Size of the title (default = 12)
#' @param x_lab Label for x-axis (default = "Comparison")
#' @param y_lab Label for y-axis (default = "Cell Type")
#' @param legend_title Title for the legend (default = "Gene Count")
#' @param x_angle Angle for x-axis text (default = 45)
#' @param top_margin Extra space at the top of the plot (default = 40)
#' @param hide_zeros Logical. If TRUE, displays zeros as NA (default = FALSE)
#' @param hide_empty_cell_types Logical. If TRUE, hides cell types that have no DEGs (default = FALSE)
#' @param min_deg_count Minimum number of DEGs required to be considered non-zero (default = 1)
#'
#' @return A ggplot object
#' @export
visualizeDEGHeatmap <- function(deg_results,
                                comparison_col = "comparison",
                                cell_type_col = "cell_type",
                                border_color = "white",
                                border_size = 0.5,
                                show_counts = TRUE,
                                count_size = 3,
                                count_color = "black",
                                low_color = "deepskyblue",
                                mid_color = "white",
                                high_color = "darkorange",
                                cluster_rows = TRUE,
                                cluster_cols = FALSE,
                                cluster_method = "complete",
                                distance_method = "euclidean",
                                group_regulation = FALSE,
                                show_reg_indicators = TRUE,
                                up_indicator_color = "darkorange",
                                down_indicator_color = "deepskyblue",
                                indicator_height = 0.2,
                                top_labels = TRUE,
                                top_label_size = 3.5,
                                label_distance = 1.5,
                                title = "DEGs Across Cell Type and Comparisons",
                                title_size = 12,
                                x_lab = "Comparison",
                                y_lab = "Cell Type",
                                legend_title = "Gene Count",
                                x_angle = 45,
                                top_margin = 40,
                                hide_zeros = FALSE,
                                hide_empty_cell_types = FALSE,
                                min_deg_count = 1) {

  # Check required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function. Please install it.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' needed for this function. Please install it.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' needed for this function. Please install it.")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' needed for this function. Please install it.")
  }

  # First, prepare the data based on input type
  if (inherits(deg_results, "list")) {
    # Extract summary data from DEG results object
    if ("summary" %in% names(deg_results)) {
      summary_data <- deg_results$summary
    } else if ("all_degs" %in% names(deg_results)) {
      # If we have all_degs, we need to create a summary
      all_degs <- deg_results$all_degs
      has_lfc_threshold_info <- "passes_lfc_threshold" %in% colnames(all_degs)

      if (has_lfc_threshold_info) {
        # Use both criteria to count DEGs
        summary_data <- all_degs %>%
          dplyr::group_by(!!dplyr::sym(cell_type_col), !!dplyr::sym(comparison_col)) %>%
          dplyr::summarize(
            n_up = sum(significant & direction == "up" & passes_lfc_threshold),
            n_down = sum(significant & direction == "down" & passes_lfc_threshold),
            .groups = "drop"
          )
      } else {
        # Fall back to just using significance
        summary_data <- all_degs %>%
          dplyr::group_by(!!dplyr::sym(cell_type_col), !!dplyr::sym(comparison_col)) %>%
          dplyr::summarize(
            n_up = sum(significant & direction == "up"),
            n_down = sum(significant & direction == "down"),
            .groups = "drop"
          )
      }
    } else {
      stop("Cannot find appropriate data in the input object")
    }
  } else if (inherits(deg_results, "data.frame")) {
    # Input is already a data frame
    has_lfc_threshold_info <- "passes_lfc_threshold" %in% colnames(deg_results)
    if (has_lfc_threshold_info) {
      summary_data <- deg_results %>%
        dplyr::group_by(!!dplyr::sym(cell_type_col), !!dplyr::sym(comparison_col)) %>%
        dplyr::summarize(
          n_up = sum(significant & direction == "up" & passes_lfc_threshold),
          n_down = sum(significant & direction == "down" & passes_lfc_threshold),
          .groups = "drop"
        )
    } else {
      summary_data <- deg_results
    }
  } else {
    stop("Input must be a data frame or a list containing DEG results")
  }

  # Ensure up and down columns exist
  if (!("n_up" %in% colnames(summary_data))) {
    if ("up_genes" %in% colnames(summary_data)) {
      summary_data$n_up <- summary_data$up_genes
    } else {
      stop("Cannot find n_up or up_genes column in data")
    }
  }

  if (!("n_down" %in% colnames(summary_data))) {
    if ("down_genes" %in% colnames(summary_data)) {
      summary_data$n_down <- summary_data$down_genes
    } else {
      stop("Cannot find n_down or down_genes column in data")
    }
  }

  # Convert small counts to zero if min_deg_count > 1
  if (min_deg_count > 1) {
    summary_data <- summary_data %>%
      dplyr::mutate(
        n_up = ifelse(n_up < min_deg_count, 0, n_up),
        n_down = ifelse(n_down < min_deg_count, 0, n_down)
      )
  }

  # Hide empty cell types if requested
  if (hide_empty_cell_types) {
    # Calculate total DEGs per cell type
    cell_type_totals <- summary_data %>%
      dplyr::group_by(!!dplyr::sym(cell_type_col)) %>%
      dplyr::summarize(total_degs = sum(n_up) + sum(n_down), .groups = "drop")

    # Keep only cell types with at least some DEGs
    non_empty_cell_types <- cell_type_totals %>%
      dplyr::filter(total_degs > 0) %>%
      dplyr::pull(!!dplyr::sym(cell_type_col))

    # Filter summary data
    summary_data <- summary_data %>%
      dplyr::filter(!!dplyr::sym(cell_type_col) %in% non_empty_cell_types)
  }

  # Check if we have data after filtering
  if (nrow(summary_data) == 0) {
    warning("No data remaining after filtering. Consider relaxing filter criteria.")

    # Create an empty plot with appropriate labels
    p <- ggplot2::ggplot() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = title,
        subtitle = "No data to display with current filtering options",
        x = x_lab,
        y = y_lab
      )

    return(p)
  }

  # Create long format data
  plot_data <- summary_data %>%
    tidyr::pivot_longer(
      cols = c("n_up", "n_down"),
      names_to = "regulation",
      values_to = "count"
    ) %>%
    dplyr::mutate(regulation = ifelse(regulation == "n_up", "Up", "Down"))

  # Apply hide_zeros if requested
  if (hide_zeros) {
    plot_data <- plot_data %>%
      dplyr::mutate(count = ifelse(count == 0, NA, count))
  }

  # Reshape data to matrix format for clustering
  wide_data <- plot_data %>%
    tidyr::unite("comparison_regulation", c(!!dplyr::sym(comparison_col), regulation), sep = "_") %>%
    tidyr::pivot_wider(
      id_cols = !!dplyr::sym(cell_type_col),
      names_from = comparison_regulation,
      values_from = count,
      values_fill = 0
    )

  # Extract row names (cell types)
  cell_types <- wide_data[[cell_type_col]]

  # Convert to matrix for clustering, excluding cell_type column
  mat_data <- as.matrix(wide_data[, -1])
  rownames(mat_data) <- cell_types

  # Replace zeros with NA for clustering if hide_zeros is TRUE
  if (hide_zeros) {
    mat_data[mat_data == 0] <- NA
  }

  # Cluster rows if requested
  if (cluster_rows && nrow(mat_data) > 1) {  # Only cluster if we have at least 2 rows
    # Handle NAs for distance calculation
    mat_for_dist <- mat_data
    mat_for_dist[is.na(mat_for_dist)] <- 0

    row_dist <- stats::dist(mat_for_dist, method = distance_method)
    row_clust <- stats::hclust(row_dist, method = cluster_method)
    cell_type_order <- rownames(mat_data)[row_clust$order]
    plot_data[[cell_type_col]] <- factor(plot_data[[cell_type_col]], levels = cell_type_order)
  } else {
    plot_data[[cell_type_col]] <- factor(plot_data[[cell_type_col]])
  }

  # Handle comparison ordering
  comparisons <- unique(plot_data[[comparison_col]])

  if (cluster_cols && length(comparisons) > 1) {  # Only cluster if we have at least 2 columns
    # Create matrix with comparison averages for clustering
    comparison_names <- comparisons
    comp_matrix <- matrix(0, nrow = length(cell_types), ncol = length(comparison_names))
    rownames(comp_matrix) <- cell_types
    colnames(comp_matrix) <- comparison_names

    for (i in 1:length(comparison_names)) {
      comp <- comparison_names[i]
      up_cols <- grepl(paste0(comp, "_Up"), colnames(mat_data))
      down_cols <- grepl(paste0(comp, "_Down"), colnames(mat_data))

      # Get column data, handling NAs
      up_data <- mat_data[, up_cols, drop = FALSE]
      down_data <- mat_data[, down_cols, drop = FALSE]

      # Replace NAs with 0 for calculation
      up_data[is.na(up_data)] <- 0
      down_data[is.na(down_data)] <- 0

      comp_matrix[, i] <- rowMeans(cbind(up_data, down_data))
    }

    col_dist <- stats::dist(t(comp_matrix), method = distance_method)
    col_clust <- stats::hclust(col_dist, method = cluster_method)
    comparison_order <- colnames(comp_matrix)[col_clust$order]
    plot_data[[comparison_col]] <- factor(plot_data[[comparison_col]], levels = comparison_order)
  } else {
    # If not clustering, still make comparison a factor to preserve order
    plot_data[[comparison_col]] <- factor(plot_data[[comparison_col]])
  }

  # Handle regulation grouping
  if (group_regulation) {
    # Create a combined x-axis position that groups all up together and all down together
    plot_data <- plot_data %>%
      dplyr::mutate(
        x_pos = ifelse(regulation == "Up",
                       as.numeric(factor(!!dplyr::sym(comparison_col))),
                       as.numeric(factor(!!dplyr::sym(comparison_col))) + length(levels(factor(!!dplyr::sym(comparison_col)))) + 1)
      )

    # Create positions for comparison labels
    up_label_positions <- 1:length(comparisons)
    down_label_positions <- (length(comparisons) + 2):(2*length(comparisons) + 1)
    all_label_positions <- c(up_label_positions, down_label_positions)

    # Create labels with line breaks for readability
    all_labels <- rep(levels(factor(plot_data[[comparison_col]])), 2)

    # Define x-axis breaks and labels
    x_breaks <- all_label_positions
    x_labels <- all_labels

    # Add indicator for Up/Down groups
    indicator_df <- data.frame(
      x = c(mean(up_label_positions), mean(down_label_positions)),
      y_min = -Inf,
      y_max = Inf,
      group = c("Up", "Down"),
      color = c(up_indicator_color, down_indicator_color)
    )

    # Create the plot
    p <- ggplot2::ggplot() +
      # Add background rectangles for up/down groups if indicators are requested
      {if(show_reg_indicators) ggplot2::geom_rect(
        data = indicator_df,
        ggplot2::aes(xmin = min(up_label_positions) - 0.5,
                     xmax = max(up_label_positions) + 0.5,
                     ymin = -Inf,
                     ymax = Inf),
        fill = up_indicator_color,
        alpha = 0.1
      )} +
      {if(show_reg_indicators) ggplot2::geom_rect(
        data = indicator_df,
        ggplot2::aes(xmin = min(down_label_positions) - 0.5,
                     xmax = max(down_label_positions) + 0.5,
                     ymin = -Inf,
                     ymax = Inf),
        fill = down_indicator_color,
        alpha = 0.1
      )} +
      # Add the main tiles
      ggplot2::geom_tile(
        data = plot_data,
        ggplot2::aes(x = x_pos, y = !!dplyr::sym(cell_type_col), fill = count),
        color = border_color,
        size = border_size,
        width = 0.9,
        height = 0.9,
        na.rm = TRUE  # Skip NA values
      ) +
      # Add count numbers if requested
      {if(show_counts) ggplot2::geom_text(
        data = plot_data,
        ggplot2::aes(x = x_pos, y = !!dplyr::sym(cell_type_col), label = count),
        size = count_size,
        color = count_color,
        na.rm = TRUE  # Skip NA values
      )} +
      # Custom color scale
      ggplot2::scale_fill_gradientn(
        colors = scales::col_numeric(
          palette = c(low_color, mid_color, high_color),
          domain = c(0, max(plot_data$count, na.rm = TRUE))
        )(seq(0, max(plot_data$count, na.rm = TRUE), length.out = 100)),
        name = legend_title,
        na.value = "grey95"  # Color for NA values
      ) +
      # Set x-axis with custom breaks and labels
      ggplot2::scale_x_continuous(
        breaks = x_breaks,
        labels = x_labels
      )

    # Add group labels at the top if requested
    if (top_labels) {
      max_y <- max(as.numeric(factor(plot_data[[cell_type_col]])))
      p <- p +
        ggplot2::annotate("text",
                          x = mean(up_label_positions),
                          y = max_y + label_distance,
                          label = "Up-regulated",
                          size = top_label_size,
                          fontface = "bold",
                          color = up_indicator_color) +
        ggplot2::annotate("text",
                          x = mean(down_label_positions),
                          y = max_y + label_distance,
                          label = "Down-regulated",
                          size = top_label_size,
                          fontface = "bold",
                          color = down_indicator_color)
    }

  } else {
    # Side-by-side approach (original)
    plot_data <- plot_data %>%
      dplyr::mutate(
        comp_id = as.numeric(factor(!!dplyr::sym(comparison_col), levels = levels(factor(!!dplyr::sym(comparison_col))))),
        x_pos = comp_id * 2 + ifelse(regulation == "Up", -0.5, 0.5)
      )

    # Create secondary x-axis labels for Up/Down
    comp_ids <- unique(plot_data$comp_id)
    up_positions <- comp_ids * 2 - 0.5
    down_positions <- comp_ids * 2 + 0.5

    # Calculate the midpoints for the up and down sections
    up_mid <- mean(range(up_positions))
    down_mid <- mean(range(down_positions))

    # Create the plot
    p <- ggplot2::ggplot()

    # Add regulation indicators at the top if requested
    if (show_reg_indicators) {
      # Define the maximum y position for consistent plotting
      max_y_level <- length(levels(factor(plot_data[[cell_type_col]])))

      # Add color indicators at the top of the plot for Up regulation
      for (pos in up_positions) {
        p <- p +
          ggplot2::annotate("rect",
                            xmin = pos - 0.45,
                            xmax = pos + 0.45,
                            ymin = max_y_level + 0.05,
                            ymax = max_y_level + 0.05 + indicator_height,
                            fill = up_indicator_color,
                            color = border_color,
                            size = border_size,
                            alpha = 0.8)
      }

      # Add color indicators at the top of the plot for Down regulation
      for (pos in down_positions) {
        p <- p +
          ggplot2::annotate("rect",
                            xmin = pos - 0.45,
                            xmax = pos + 0.45,
                            ymin = max_y_level + 0.05,
                            ymax = max_y_level + 0.05 + indicator_height,
                            fill = down_indicator_color,
                            color = border_color,
                            size = border_size,
                            alpha = 0.8)
      }

      # Add Up/Down labels at the top if requested, positioned further above
      if (top_labels) {
        # Add labels for entire Up/Down regions, centered
        p <- p +
          ggplot2::annotate("text",
                            x = up_mid,
                            y = max_y_level + label_distance,
                            label = "Up-regulated",
                            size = top_label_size,
                            fontface = "bold",
                            color = up_indicator_color) +
          ggplot2::annotate("text",
                            x = down_mid,
                            y = max_y_level + label_distance,
                            label = "Down-regulated",
                            size = top_label_size,
                            fontface = "bold",
                            color = down_indicator_color)
      }
    }

    # Add dummy points with regulation colors for the legend
    reg_colors <- c("Up" = up_indicator_color, "Down" = down_indicator_color)
    dummy_data <- data.frame(
      x = rep(0, 2),
      y = rep(0, 2),
      regulation = c("Up", "Down")
    )

    p <- p +
      # Add invisible points for regulation legend
      ggplot2::geom_point(
        data = dummy_data,
        ggplot2::aes(x = x, y = y, color = regulation),
        shape = 15,   # Square shape for better legend appearance
        size = 5,     # Larger size for visibility in legend
        alpha = 0     # Make invisible in the plot
      ) +
      # Define regulation colors
      ggplot2::scale_color_manual(
        values = reg_colors,
        name = "Regulation"
      )

    # Add the main tiles
    p <- p +
      ggplot2::geom_tile(
        data = plot_data,
        ggplot2::aes(x = x_pos, y = !!dplyr::sym(cell_type_col), fill = count),
        color = border_color,
        size = border_size,
        width = 0.9,
        height = 0.9,
        na.rm = TRUE  # Skip NA values
      )

    # Add count numbers if requested
    if (show_counts) {
      p <- p + ggplot2::geom_text(
        data = plot_data,
        ggplot2::aes(x = x_pos, y = !!dplyr::sym(cell_type_col), label = count),
        size = count_size,
        color = count_color,
        na.rm = TRUE  # Skip NA values
      )
    }

    # Custom color scale
    p <- p + ggplot2::scale_fill_gradientn(
      colors = scales::col_numeric(
        palette = c(low_color, mid_color, high_color),
        domain = c(0, max(plot_data$count, na.rm = TRUE))
      )(seq(0, max(plot_data$count, na.rm = TRUE), length.out = 100)),
      name = legend_title,
      na.value = "grey95"  # Color for NA values
    )

    # Set x-axis with labels for comparisons (no secondary axis)
    p <- p + ggplot2::scale_x_continuous(
      breaks = comp_ids * 2,
      labels = levels(factor(plot_data[[comparison_col]]))
    )
  }

  # Define the maximum y position for expansion, based on label distance
  max_y_expansion <- if(top_labels) label_distance + 0.5 else indicator_height + 0.3

  # Complete the plot styling
  p <- p +
    ggplot2::labs(
      title = title,
      x = x_lab,
      y = y_lab
    ) +
    # Expand y-axis to make room for the top labels
    ggplot2::scale_y_discrete(
      expand = ggplot2::expansion(add = c(0, max_y_expansion))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_size, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = x_angle, hjust = 1),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = top_margin, r = 20, b = 20, l = 20, unit = "pt"),
      legend.position = "right"
    )

  return(p)
}

#' Create a Barplot of Gene Regulation by Cell Type
#'
#' Creates a barplot showing up and down regulated genes for each cell type
#' across different comparisons.
#'
#' @param deg_results Output from findDifferentialGenes or compareDEGs function
#' @param comparison_col Column name for the comparison categories
#' @param cell_type_col Column name for the cell types
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
#' @param x_lab Label for x-axis (default = NULL, will be automatically set)
#' @param y_lab Label for y-axis (default = NULL, will be automatically set)
#' @param legend_title Title for the legend (default = "Regulation")
#' @param x_angle Angle for x-axis text (default = 45)
#' @param hide_zeros Logical. If TRUE, hides bars for cell types with zero DEGs (default = FALSE)
#' @param hide_empty_cell_types Logical. If TRUE, hides cell types that have no DEGs in any comparison (default = FALSE)
#' @param min_deg_count Minimum number of DEGs required to show a bar (default = 0)
#'
#' @return A ggplot object
#' @export
visualizeDEGBarplot <- function(deg_results,
                                comparison_col = "comparison",
                                cell_type_col = "cell_type",
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
                                x_angle = 45,
                                hide_zeros = FALSE,
                                hide_empty_cell_types = FALSE,
                                min_deg_count = 0) {

  # Check required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function. Please install it.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' needed for this function. Please install it.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' needed for this function. Please install it.")
  }

  # First, prepare the data based on input type
  if (inherits(deg_results, "list")) {
    # Extract summary data from DEG results object
    if ("summary" %in% names(deg_results)) {
      summary_data <- deg_results$summary
    } else if ("all_degs" %in% names(deg_results)) {
      # If we have all_degs, we need to create a summary
      all_degs <- deg_results$all_degs

      has_lfc_threshold_info <- "passes_lfc_threshold" %in% colnames(all_degs)

      if (has_lfc_threshold_info) {
        # Use both criteria to count DEGs
        summary_data <- all_degs %>%
          dplyr::group_by(!!dplyr::sym(cell_type_col), !!dplyr::sym(comparison_col)) %>%
          dplyr::summarize(
            n_up = sum(significant & direction == "up" & passes_lfc_threshold),
            n_down = sum(significant & direction == "down" & passes_lfc_threshold),
            .groups = "drop"
          )
      } else {
        # Fall back to just using significance
        summary_data <- all_degs %>%
          dplyr::group_by(!!dplyr::sym(cell_type_col), !!dplyr::sym(comparison_col)) %>%
          dplyr::summarize(
            n_up = sum(significant & direction == "up"),
            n_down = sum(significant & direction == "down"),
            .groups = "drop"
          )
      }
    } else {
      stop("Cannot find appropriate data in the input object")
    }
  } else if (inherits(deg_results, "data.frame")) {
    # Input is already a data frame
    has_lfc_threshold_info <- "passes_lfc_threshold" %in% colnames(deg_results)

    if (has_lfc_threshold_info) {
      summary_data <- deg_results %>%
        dplyr::group_by(!!dplyr::sym(cell_type_col), !!dplyr::sym(comparison_col)) %>%
        dplyr::summarize(
          n_up = sum(significant & direction == "up" & passes_lfc_threshold),
          n_down = sum(significant & direction == "down" & passes_lfc_threshold),
          .groups = "drop"
        )
    } else {
      summary_data <- deg_results
    }
    } else {
    stop("Input must be a data frame or a list containing DEG results")
  }

  # Ensure up and down columns exist
  if (!("n_up" %in% colnames(summary_data))) {
    if ("up_genes" %in% colnames(summary_data)) {
      summary_data$n_up <- summary_data$up_genes
    } else {
      stop("Cannot find n_up or up_genes column in data")
    }
  }

  if (!("n_down" %in% colnames(summary_data))) {
    if ("down_genes" %in% colnames(summary_data)) {
      summary_data$n_down <- summary_data$down_genes
    } else {
      stop("Cannot find n_down or down_genes column in data")
    }
  }

  # Validate facet_by parameter
  if (!facet_by %in% c("comparison", "cell_type", "none")) {
    stop("facet_by must be one of: 'comparison', 'cell_type', or 'none'")
  }

  # Hide empty cell types if requested
  if (hide_empty_cell_types) {
    # Calculate total DEGs per cell type
    cell_type_totals <- summary_data %>%
      dplyr::group_by(!!dplyr::sym(cell_type_col)) %>%
      dplyr::summarize(total_degs = sum(n_up) + sum(n_down), .groups = "drop")

    # Keep only cell types with at least some DEGs
    non_empty_cell_types <- cell_type_totals %>%
      dplyr::filter(total_degs > 0) %>%
      dplyr::pull(!!dplyr::sym(cell_type_col))

    # Filter summary data
    summary_data <- summary_data %>%
      dplyr::filter(!!dplyr::sym(cell_type_col) %in% non_empty_cell_types)
  }

  # Reshape data to long format for plotting
  plot_data_long <- summary_data %>%
    tidyr::pivot_longer(
      cols = c("n_up", "n_down"),
      names_to = "regulation",
      values_to = "count"
    ) %>%
    dplyr::mutate(
      regulation = ifelse(regulation == "n_up", "Up", "Down"),
      # Make sure cell_type and comparison_pattern are factors to preserve order
      !!dplyr::sym(cell_type_col) := factor(!!dplyr::sym(cell_type_col)),
      !!dplyr::sym(comparison_col) := factor(!!dplyr::sym(comparison_col))
    )

  # Filter out zeros if requested
  if (hide_zeros) {
    plot_data_long <- plot_data_long %>%
      dplyr::filter(count > min_deg_count)
  }

  # Check if we have data after filtering
  if (nrow(plot_data_long) == 0) {
    warning("No data remaining after filtering. Consider relaxing filter criteria.")

    # Create an empty plot with appropriate labels
    p <- ggplot2::ggplot() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = title,
        subtitle = "No data to display with current filtering options",
        x = if (is.null(x_lab)) (if (horizontal) "Gene Count" else "Cell Type") else x_lab,
        y = if (is.null(y_lab)) (if (horizontal) "Cell Type" else "Gene Count") else y_lab
      )

    return(p)
  }

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
    p <- ggplot2::ggplot(plot_data_long) +
      ggplot2::geom_col(
        ggplot2::aes(
          y = !!dplyr::sym(cell_type_col),
          x = count,
          fill = regulation
        ),
        position = position_type,
        width = bar_width
      )
  } else {
    # Vertical bars
    p <- ggplot2::ggplot(plot_data_long) +
      ggplot2::geom_col(
        ggplot2::aes(
          x = !!dplyr::sym(cell_type_col),
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
        p <- p + ggplot2::geom_text(
          ggplot2::aes(
            y = !!dplyr::sym(cell_type_col),
            x = count/2,
            label = count,
            group = regulation
          ),
          position = ggplot2::position_stack(),
          size = count_size,
          color = count_color
        )
      } else {
        # For dodged bars, position text at the end of each bar
        p <- p + ggplot2::geom_text(
          ggplot2::aes(
            y = !!dplyr::sym(cell_type_col),
            x = count,
            label = count,
            group = regulation
          ),
          position = ggplot2::position_dodge(width = bar_width),
          hjust = -0.2,
          size = count_size,
          color = count_color
        )
      }
    } else {
      # For vertical bars
      if (stacked) {
        # For stacked bars, position text in the middle of each segment
        p <- p + ggplot2::geom_text(
          ggplot2::aes(
            x = !!dplyr::sym(cell_type_col),
            y = count/2,
            label = count,
            group = regulation
          ),
          position = ggplot2::position_stack(),
          size = count_size,
          color = count_color
        )
      } else {
        # For dodged bars, position text at the top of each bar
        p <- p + ggplot2::geom_text(
          ggplot2::aes(
            x = !!dplyr::sym(cell_type_col),
            y = count,
            label = count,
            group = regulation
          ),
          position = ggplot2::position_dodge(width = bar_width),
          vjust = -0.2,
          size = count_size,
          color = count_color
        )
      }
    }
  }

  # Add fill colors and legend
  p <- p + ggplot2::scale_fill_manual(
    values = c("Up" = up_color, "Down" = down_color),
    name = legend_title
  )

  # Add facets if requested
  if (facet_by != "none") {
    if (facet_by == "comparison") {
      p <- p + ggplot2::facet_wrap(
        ggplot2::vars(!!dplyr::sym(comparison_col)),
        scales = facet_scales,
        ncol = facet_ncol
      )
    } else if (facet_by == "cell_type") {
      p <- p + ggplot2::facet_wrap(
        ggplot2::vars(!!dplyr::sym(cell_type_col)),
        scales = facet_scales,
        ncol = facet_ncol
      )
    }
  }

  # Complete the plot styling
  p <- p + ggplot2::labs(
    title = title,
    x = x_lab,
    y = y_lab
  ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_size, face = "bold"),
      axis.text.x = if (!horizontal) ggplot2::element_text(angle = x_angle, hjust = 1) else ggplot2::element_text(),
      legend.position = "right",
      panel.grid.major.x = if (horizontal) ggplot2::element_line(color = "grey90") else ggplot2::element_blank(),
      panel.grid.major.y = if (!horizontal) ggplot2::element_line(color = "grey90") else ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.text = ggplot2::element_text(face = "bold")
    )

  return(p)
}

#' Create a Gene Regulation Dotplot
#'
#' Creates a dotplot visualizing up and down regulated genes across different comparisons
#' with customizable shapes and count label options.
#'
#' @param deg_results Output from findDifferentialGenes or compareDEGs function
#' @param comparison_col Column name for the comparison categories
#' @param cell_type_col Column name for the cell types
#' @param count_display Where to display count numbers: "none", "inside", "below", or "both" (default = "none")
#' @param use_triangles Logical. If TRUE, uses triangles (shapes 24, 25), otherwise uses circles (shape 21) (default = TRUE)
#' @param point_size_range Numeric vector of length 2 for the range of point sizes (default = c(1, 10))
#' @param up_color Color for up-regulated genes (default = "deepskyblue")
#' @param down_color Color for down-regulated genes (default = "darkorange")
#' @param offset Size of horizontal offset between up/down points (default = 0.2)
#' @param title Plot title (default = "DEGs Across Cell Type and Comparisons")
#' @param x_angle Angle for x-axis text (default = 45)
#' @param hide_zeros Logical. If TRUE, hides dots for cell types with zero DEGs (default = FALSE)
#' @param hide_empty_cell_types Logical. If TRUE, hides cell types that have no DEGs in any comparison (default = FALSE)
#' @param min_deg_count Minimum number of DEGs required to show a dot (default = 0)
#'
#' @return A ggplot object
#' @export
visualizeDEGDotplot <- function(deg_results,
                                comparison_col = "comparison",
                                cell_type_col = "cell_type",
                                count_display = "none",
                                use_triangles = TRUE,
                                point_size_range = c(1, 10),
                                up_color = "deepskyblue",
                                down_color = "darkorange",
                                offset = 0.2,
                                title = "DEGs Across Cell Type and Comparisons",
                                x_angle = 45,
                                hide_zeros = FALSE,
                                hide_empty_cell_types = FALSE,
                                min_deg_count = 0) {

  # Check required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function. Please install it.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' needed for this function. Please install it.")
  }

  # Validate count_display parameter
  if (!count_display %in% c("none", "inside", "below", "both")) {
    stop("count_display must be one of: 'none', 'inside', 'below', or 'both'")
  }

  # First, prepare the data based on input type
  if (inherits(deg_results, "list")) {
    # Extract summary data from DEG results object
    if ("summary" %in% names(deg_results)) {
      summary_data <- deg_results$summary
    } else if ("all_degs" %in% names(deg_results)) {
      # If we have all_degs, we need to create a summary
      has_lfc_threshold_info <- "passes_lfc_threshold" %in% colnames(all_degs)

      if (has_lfc_threshold_info) {
        # Create summary with both criteria
        plot_data <- all_degs %>%
          dplyr::filter(significant & passes_lfc_threshold) %>%
          dplyr::group_by(!!dplyr::sym(cell_type_col), !!dplyr::sym(comparison_col), direction) %>%
          dplyr::summarize(count = n(), .groups = "drop") %>%
          tidyr::pivot_wider(
            names_from = direction,
            values_from = count,
            values_fill = 0,
            names_prefix = "n_"
          )

        # Ensure we have both up and down columns
        if (!"n_up" %in% names(plot_data)) plot_data$n_up <- 0
        if (!"n_down" %in% names(plot_data)) plot_data$n_down <- 0

        summary_data <- plot_data
      } else {
        # Just use significance
        summary_data <- all_degs %>%
          dplyr::filter(significant) %>%
          dplyr::group_by(!!dplyr::sym(cell_type_col), !!dplyr::sym(comparison_col), direction) %>%
          dplyr::summarize(count = n(), .groups = "drop") %>%
          tidyr::pivot_wider(
            names_from = direction,
            values_from = count,
            values_fill = 0,
            names_prefix = "n_"
          )

        # Ensure we have both up and down columns
        if (!"n_up" %in% names(summary_data)) summary_data$n_up <- 0
        if (!"n_down" %in% names(summary_data)) summary_data$n_down <- 0
      }
    } else {
      stop("Cannot find appropriate data in the input object")
    }
  } else if (inherits(deg_results, "data.frame")) {
    # Input is already a data frame
    has_lfc_threshold_info <- "passes_lfc_threshold" %in% colnames(deg_results)

    if (has_lfc_threshold_info) {
      summary_data <- deg_results %>%
        dplyr::group_by(!!dplyr::sym(cell_type_col), !!dplyr::sym(comparison_col)) %>%
        dplyr::summarize(
          n_up = sum(significant & direction == "up" & passes_lfc_threshold),
          n_down = sum(significant & direction == "down" & passes_lfc_threshold),
          .groups = "drop"
        )
    } else {
      summary_data <- deg_results
    }
    } else {
    stop("Input must be a data frame or a list containing DEG results")
  }

  # Ensure up and down columns exist
  if (!("n_up" %in% colnames(summary_data))) {
    if ("up_genes" %in% colnames(summary_data)) {
      summary_data$n_up <- summary_data$up_genes
    } else {
      stop("Cannot find n_up or up_genes column in data")
    }
  }

  if (!("n_down" %in% colnames(summary_data))) {
    if ("down_genes" %in% colnames(summary_data)) {
      summary_data$n_down <- summary_data$down_genes
    } else {
      stop("Cannot find n_down or down_genes column in data")
    }
  }

  # Reshape data to long format
  plot_data <- summary_data %>%
    tidyr::pivot_longer(
      cols = c("n_up", "n_down"),
      names_to = "regulation",
      values_to = "count"
    ) %>%
    dplyr::mutate(regulation = ifelse(regulation == "n_up", "Up", "Down"))

  # Filter out zeros if requested
  if (hide_zeros) {
    plot_data <- plot_data %>%
      dplyr::filter(count > min_deg_count)
  }

  # Hide empty cell types if requested
  if (hide_empty_cell_types) {
    # Calculate total DEGs per cell type
    cell_type_totals <- plot_data %>%
      dplyr::group_by(!!dplyr::sym(cell_type_col)) %>%
      dplyr::summarize(total_degs = sum(count), .groups = "drop")

    # Keep only cell types with at least some DEGs
    non_empty_cell_types <- cell_type_totals %>%
      dplyr::filter(total_degs > 0) %>%
      dplyr::pull(!!dplyr::sym(cell_type_col))

    # Filter plot data
    plot_data <- plot_data %>%
      dplyr::filter(!!dplyr::sym(cell_type_col) %in% non_empty_cell_types)
  }

  # Check if we have data after filtering
  if (nrow(plot_data) == 0) {
    warning("No data remaining after filtering. Consider relaxing filter criteria.")

    # Create an empty plot with appropriate labels
    p <- ggplot2::ggplot() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = title,
        subtitle = "No data to display with current filtering options",
        x = "Comparison",
        y = "Cell Type"
      )

    return(p)
  }

  # Create numeric position variable for the x-axis
  plot_data <- plot_data %>%
    # Get unique comparisons
    dplyr::mutate(comp_id = as.numeric(factor(!!dplyr::sym(comparison_col)))) %>%
    # Add offset for regulation (Up/Down)
    dplyr::mutate(x_pos = comp_id + ifelse(regulation == "Up", -offset, offset)) %>%
    # Convert cell_type to a factor
    dplyr::mutate(!!dplyr::sym(cell_type_col) := factor(!!dplyr::sym(cell_type_col)))

  # Set shapes based on the use_triangles parameter
  if (use_triangles) {
    # Triangle shapes (up/down)
    shape_values <- c("Up" = 24, "Down" = 25)

    # Create the base plot with triangles
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(
      x = x_pos,
      y = !!dplyr::sym(cell_type_col),
      size = count,
      color = regulation,
      shape = regulation
    )) +
      ggplot2::geom_point(alpha = 0.9, stroke = 2, fill = "white") +
      ggplot2::scale_shape_manual(values = shape_values)
  } else {
    # Circle shape for both, color differentiates up/down
    # Create the base plot with circles
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(
      x = x_pos,
      y = !!dplyr::sym(cell_type_col),
      size = count,
      color = regulation
    )) +
      ggplot2::geom_point(alpha = 0.9, shape = 21, stroke = 2) +
      ggplot2::scale_color_manual(values = c("Up" = ggplot2::alpha(up_color, 0.3), "Down" = ggplot2::alpha(down_color, 0.3)))
  }

  # Add common aesthetic elements
  p <- p +
    ggplot2::scale_color_manual(values = c("Up" = up_color, "Down" = down_color)) +
    ggplot2::scale_size_continuous(range = point_size_range) +
    # Use unique comparison names for x-axis labels
    ggplot2::scale_x_continuous(
      breaks = unique(plot_data$comp_id),
      labels = unique(plot_data[[comparison_col]])
    ) +
    ggplot2::labs(
      title = title,
      x = "Comparison",
      y = "Cell Type",
      size = "Gene Count",
      color = "Regulation"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = x_angle, hjust = 1),
      panel.grid.major = ggplot2::element_line(color = "grey90"),
      panel.grid.minor = ggplot2::element_line(color = "grey95"),
      legend.position = "right"
    )

  # Add shape to legend if using triangles
  if (use_triangles) {
    p <- p + ggplot2::labs(shape = "Regulation")
  } else {
    p <- p + ggplot2::labs(fill = "Regulation")
  }

  # Add text based on count_display parameter
  if (count_display %in% c("inside", "both")) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = count),
                                size = 2.5,
                                show.legend = FALSE,
                                color = "black")
  }

  if (count_display %in% c("below", "both")) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = count),
                                nudge_y = -0.3,
                                size = 3,
                                show.legend = FALSE)
  }

  return(p)
}


#' Visualize DEG Patterns Across Cell Types and Conditions
#'
#' Creates a multi-panel visualization combining different plot types for comprehensive
#' gene regulation pattern analysis.
#'
#' @param deg_results Output from findDifferentialGenes or compareDEGs function
#' @param plot_types Vector of plot types to include: "heatmap", "barplot", "dotplot"
#' @param ncol Number of columns for the multi-panel layout (default = 1)
#' @param comparison_col Column name for the comparison categories
#' @param cell_type_col Column name for the cell types
#' @param up_color Color for up-regulated genes (default = "darkorange")
#' @param down_color Color for down-regulated genes (default = "deepskyblue")
#' @param title Overall title for the combined plot
#' @param filter_min_degs Minimum number of DEGs required to include a cell type-comparison
#' @param cluster_rows Logical. Whether to cluster rows in heatmap
#' @param hide_zeros Logical. If TRUE, hides dots/bars for cell types with zero DEGs
#' @param hide_empty_cell_types Logical. If TRUE, hides cell types that have no DEGs
#'
#' @return A combined plot (patchwork object)
#' @export
visualizeDEGPatterns <- function(deg_results,
                                 plot_types = c("heatmap", "barplot"),
                                 ncol = 1,
                                 comparison_col = "comparison",
                                 cell_type_col = "cell_type",
                                 up_color = "darkorange",
                                 down_color = "deepskyblue",
                                 title = "DEG Patterns Across Cell Types and Conditions",
                                 filter_min_degs = 0,
                                 cluster_rows = TRUE,
                                 hide_zeros = FALSE,
                                 hide_empty_cell_types = FALSE) {

  # Check required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function. Please install it.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' needed for this function. Please install it.")
  }

  # Validate plot_types
  valid_types <- c("heatmap", "barplot", "dotplot")
  if (!all(plot_types %in% valid_types)) {
    stop("All plot_types must be one of: 'heatmap', 'barplot', or 'dotplot'")
  }

  # First, prepare the data based on input type
  if (inherits(deg_results, "list")) {
    # Extract summary data from DEG results object
    if ("summary" %in% names(deg_results)) {
      summary_data <- deg_results$summary
    } else if ("all_degs" %in% names(deg_results)) {
      # If we have all_degs, we need to create a summary
      all_degs <- deg_results$all_degs

      # Check if we have the passes_lfc_threshold column
      has_lfc_threshold_info <- "passes_lfc_threshold" %in% colnames(all_degs)

      if (has_lfc_threshold_info) {
        # Generate a modified deg_results with the correct summary
        summary_data <- all_degs %>%
          dplyr::group_by(!!dplyr::sym(cell_type_col), !!dplyr::sym(comparison_col)) %>%
          dplyr::summarize(
            n_up = sum(significant & direction == "up" & passes_lfc_threshold),
            n_down = sum(significant & direction == "down" & passes_lfc_threshold),
            n_genes_tested = n(),
            n_significant = sum(significant & passes_lfc_threshold),
            .groups = "drop"
          )
      } else {
        summary_data <- all_degs %>%
          dplyr::group_by(!!dplyr::sym(cell_type_col), !!dplyr::sym(comparison_col)) %>%
          dplyr::summarize(
            n_up = sum(significant & direction == "up"),
            n_down = sum(significant & direction == "down"),
            n_genes_tested = n(),
            n_significant = sum(significant),
            .groups = "drop"
          )
      }
    } else {
      stop("Cannot find appropriate data in the input object")
    }
  } else if (inherits(deg_results, "data.frame")) {
    # Input is already a data frame
    summary_data <- deg_results
  } else {
    stop("Input must be a data frame or a list containing DEG results")
  }

  # Ensure up and down columns exist
  if (!("n_up" %in% colnames(summary_data))) {
    if ("up_genes" %in% colnames(summary_data)) {
      summary_data$n_up <- summary_data$up_genes
    } else {
      stop("Cannot find n_up or up_genes column in data")
    }
  }

  if (!("n_down" %in% colnames(summary_data))) {
    if ("down_genes" %in% colnames(summary_data)) {
      summary_data$n_down <- summary_data$down_genes
    } else {
      stop("Cannot find n_down or down_genes column in data")
    }
  }

  # Convert small counts to zero if filter_min_degs > 0
  if (filter_min_degs > 0) {
    summary_data <- summary_data %>%
      dplyr::mutate(
        n_up = ifelse(n_up < filter_min_degs, 0, n_up),
        n_down = ifelse(n_down < filter_min_degs, 0, n_down)
      )
  }

  # Hide empty cell types if requested
  if (hide_empty_cell_types) {
    # Calculate total DEGs per cell type
    cell_type_totals <- summary_data %>%
      dplyr::group_by(!!dplyr::sym(cell_type_col)) %>%
      dplyr::summarize(total_degs = sum(n_up) + sum(n_down), .groups = "drop")

    # Keep only cell types with at least some DEGs
    non_empty_cell_types <- cell_type_totals %>%
      dplyr::filter(total_degs > 0) %>%
      dplyr::pull(!!dplyr::sym(cell_type_col))

    # Filter summary data
    summary_data <- summary_data %>%
      dplyr::filter(!!dplyr::sym(cell_type_col) %in% non_empty_cell_types)
  }

  # Check if we have data after filtering
  if (nrow(summary_data) == 0) {
    warning("No data remaining after filtering. Consider relaxing filter criteria.")

    # Create an empty plot with appropriate labels
    p <- ggplot2::ggplot() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = title,
        subtitle = "No data to display with current filtering options",
        x = "Comparison",
        y = "Cell Type"
      )

    return(p)
  }

  # Create individual plots
  plots <- list()

  for (plot_type in plot_types) {
    if (plot_type == "heatmap") {
      plots[["heatmap"]] <- visualizeDEGHeatmap(
        deg_results = summary_data,
        comparison_col = comparison_col,
        cell_type_col = cell_type_col,
        up_indicator_color = up_color,
        down_indicator_color = down_color,
        cluster_rows = cluster_rows,
        title = "DEG Heatmap"
      )
    } else if (plot_type == "barplot") {
      # Create long format data for barplot with filtering
      barplot_data <- summary_data

      # Filter out zeros if requested (for barplot)
      if (hide_zeros) {
        barplot_data <- barplot_data %>%
          dplyr::mutate(
            n_up = ifelse(n_up == 0, NA, n_up),
            n_down = ifelse(n_down == 0, NA, n_down)
          )
      }

      plots[["barplot"]] <- visualizeDEGBarplot(
        deg_results = barplot_data,
        comparison_col = comparison_col,
        cell_type_col = cell_type_col,
        up_color = up_color,
        down_color = down_color,
        facet_by = "comparison",
        title = "DEG Counts by Cell Type"
      )
    } else if (plot_type == "dotplot") {
      # Create dotplot
      plots[["dotplot"]] <- visualizeDEGDotplot(
        deg_results = summary_data,
        comparison_col = comparison_col,
        cell_type_col = cell_type_col,
        up_color = up_color,
        down_color = down_color,
        count_display = "inside",
        title = "DEG Dotplot",
        hide_zeros = hide_zeros,
        hide_empty_cell_types = hide_empty_cell_types
      )
    }
  }

  # Combine plots using patchwork
  combined_plot <- patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5)
      )
    )

  return(combined_plot)
}
