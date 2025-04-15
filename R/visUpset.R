# upset.R

#' Create an UpSet Plot for Visualizing Set Intersections
#'
#' Creates a custom UpSet plot showing the intersections between sets,
#' particularly useful for visualizing shared DEGs across cell types or comparisons.
#'
#' @param data_list A named list of vectors, each containing elements in a set
#' @param sets Optional vector of set names to include (default: all sets in data_list)
#' @param min_intersection_size Minimum intersection size to include (default: 1)
#' @param max_sets_display Maximum number of sets to display (default: all)
#' @param sort_sets_by How to sort the sets: "size", "name", or "custom" (default: "size")
#' @param sort_sets_decreasing Whether to sort sets in decreasing order (default: TRUE)
#' @param custom_sets_order Custom order for sets if sort_sets_by="custom"
#' @param sort_intersections_by How to sort intersections: "freq", "degree", or "custom" (default: "freq")
#' @param sort_intersections_decreasing Whether to sort intersections in decreasing order (default: TRUE)
#' @param custom_intersections_order Custom order for intersections if sort_intersections_by="custom"
#' @param intersection_color Color for intersection dots and lines (default: "black")
#' @param main_bar_color Color for the intersection size bars (default: "steelblue")
#' @param sets_bar_colors Named vector of colors for each set (default: auto-generated)
#' @param highlight_intersections Vector of intersection IDs to highlight (default: NULL)
#' @param highlight_color Color for highlighted intersections (default: "darkorange")
#' @param point_size Size of points in the matrix (default: 3)
#' @param line_size Width of connecting lines (default: 1)
#' @param bar_width Width of bars (0-1 scale) (default: 0.7)
#' @param text_angle Angle for text labels (default: 0)
#' @param text_size Size of text in the plot (default: 10)
#' @param title Plot title (default: NULL)
#' @param show_numbers Whether to show numbers on bars (default: TRUE)
#' @param intersection_title Title for the intersection size plot (default: "Intersection Size")
#' @param set_size_title Title for the set size plot (default: "Set Size")
#'
#' @return A ggplot2 object or a combined patchwork layout
#' @export
upsetPlot <- function(data_list,
                      sets = NULL,
                      min_intersection_size = 1,
                      max_sets_display = NULL,
                      sort_sets_by = "size",
                      sort_sets_decreasing = TRUE,
                      custom_sets_order = NULL,
                      sort_intersections_by = "freq",
                      sort_intersections_decreasing = TRUE,
                      custom_intersections_order = NULL,
                      intersection_color = "black",
                      main_bar_color = "steelblue",
                      sets_bar_colors = NULL,
                      highlight_intersections = NULL,
                      highlight_color = "darkorange",
                      point_size = 3,
                      line_size = 1,
                      bar_width = 0.7,
                      text_angle = 0,
                      text_size = 10,
                      title = NULL,
                      show_numbers = TRUE,
                      intersection_title = "Intersection Size",
                      set_size_title = "Set Size") {

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
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' needed for this function. Please install it.")
  }

  # Parameter validation
  if (!is.list(data_list)) {
    stop("'data_list' must be a list where each element is a vector of elements in that set")
  }

  # If sets parameter is not provided, use all names from the list
  if (is.null(sets)) {
    sets <- names(data_list)
    if (is.null(sets)) {
      sets <- paste0("Set", 1:length(data_list))
      names(data_list) <- sets
    }
  } else {
    # Validate if all specified sets exist in data_list
    if (!all(sets %in% names(data_list))) {
      stop("Not all specified sets exist in data_list")
    }
    # Filter data_list to only include specified sets
    data_list <- data_list[sets]
  }

  # Limit the number of sets to display if specified
  if (!is.null(max_sets_display) && max_sets_display < length(sets)) {
    set_sizes <- sapply(data_list, length)
    sets_to_keep <- if (sort_sets_by == "size") {
      names(sort(set_sizes, decreasing = sort_sets_decreasing)[1:max_sets_display])
    } else if (sort_sets_by == "name") {
      sort(names(set_sizes), decreasing = sort_sets_decreasing)[1:max_sets_display]
    } else {
      sets[1:max_sets_display]
    }

    data_list <- data_list[sets_to_keep]
    sets <- sets_to_keep
  }

  # Get all unique elements
  all_elements <- unique(unlist(data_list))

  # Create a binary membership matrix
  membership_df <- data.frame(element = all_elements)

  for (set_name in sets) {
    membership_df[[set_name]] <- as.integer(membership_df$element %in% data_list[[set_name]])
  }

  # Generate all possible intersections by grouping and counting
  pattern_df <- membership_df %>%
    dplyr::select(-element) %>%
    dplyr::group_by_all() %>%
    dplyr::summarise(size = dplyr::n(), .groups = "drop")

  # Filter by minimum intersection size
  pattern_df <- pattern_df %>%
    dplyr::filter(size >= min_intersection_size)

  if (nrow(pattern_df) == 0) {
    stop("No intersections meet the minimum size requirement")
  }

  # Add degree (number of sets in each intersection)
  pattern_df$degree <- rowSums(dplyr::select(pattern_df, -size))

  # Sort intersections
  if (sort_intersections_by == "freq") {
    pattern_df <- pattern_df %>%
      dplyr::arrange(if(sort_intersections_decreasing) dplyr::desc(size) else size)
  } else if (sort_intersections_by == "degree") {
    pattern_df <- pattern_df %>%
      dplyr::arrange(
        if(sort_intersections_decreasing) dplyr::desc(degree) else degree,
        if(sort_intersections_decreasing) dplyr::desc(size) else size
      )
  } else if (sort_intersections_by == "custom" && !is.null(custom_intersections_order)) {
    # Use custom order - this assumes custom_intersections_order contains pattern IDs in desired order
    pattern_df <- pattern_df[custom_intersections_order, ]
  }

  # Create a pattern ID for each intersection
  pattern_df$pattern_id <- 1:nrow(pattern_df)

  # Calculate set sizes
  set_sizes <- sapply(data_list, length)
  set_size_df <- data.frame(
    set_name = names(set_sizes),
    size = as.numeric(set_sizes)
  )

  # Sort sets based on specified order
  if (sort_sets_by == "size") {
    sets_ordered <- set_size_df %>%
      dplyr::arrange(if(sort_sets_decreasing) dplyr::desc(size) else size) %>%
      dplyr::pull(set_name)

    set_size_df$set_name <- factor(set_size_df$set_name, levels = sets_ordered)
    sets <- sets_ordered
  } else if (sort_sets_by == "name") {
    sets_ordered <- sort(sets, decreasing = sort_sets_decreasing)
    set_size_df$set_name <- factor(set_size_df$set_name, levels = sets_ordered)
    sets <- sets_ordered
  } else if (sort_sets_by == "custom" && !is.null(custom_sets_order)) {
    # Verify custom order contains all sets
    if (!all(sets %in% custom_sets_order)) {
      warning("custom_sets_order doesn't contain all sets. Using size ordering instead.")
      sets_ordered <- set_size_df %>%
        dplyr::arrange(dplyr::desc(size)) %>%
        dplyr::pull(set_name)
    } else {
      sets_ordered <- custom_sets_order[custom_sets_order %in% sets]
    }
    set_size_df$set_name <- factor(set_size_df$set_name, levels = sets_ordered)
    sets <- sets_ordered
  }

  # Convert set names to numeric positions for uniform positioning
  set_size_df$y_pos <- 1:nrow(set_size_df)

  # Convert pattern_df to long format for matrix visualization
  pattern_long <- pattern_df %>%
    tidyr::pivot_longer(cols = dplyr::all_of(sets),
                        names_to = "set_name",
                        values_to = "is_in_set") %>%
    dplyr::select(pattern_id, set_name, is_in_set, size, degree)

  # Ensure consistent factor levels for sets
  pattern_long$set_name <- factor(pattern_long$set_name, levels = sets)

  # Join with set positions
  pattern_long <- pattern_long %>%
    dplyr::left_join(set_size_df %>% dplyr::select(set_name, y_pos), by = "set_name")

  # Define set colors if not provided
  if (is.null(sets_bar_colors)) {
    # Default color palette
    default_colors <- c(
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
      "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
      "#999999", "#66C2A5", "#FC8D62", "#8DA0CB"
    )

    # Generate colors for all sets
    sets_bar_colors <- stats::setNames(
      default_colors[1:min(length(sets), length(default_colors))],
      sets
    )

    # For any remaining sets beyond the default colors, recycle colors
    if (length(sets) > length(default_colors)) {
      remaining_sets <- sets[(length(default_colors) + 1):length(sets)]
      remaining_colors <- default_colors[1:length(remaining_sets)]
      additional_colors <- stats::setNames(remaining_colors, remaining_sets)
      sets_bar_colors <- c(sets_bar_colors, additional_colors)
    }
  }

  # Create a position for each intersection
  pattern_df$x_pos <- 1:nrow(pattern_df)

  # Add x_pos to pattern_long for matrix plot
  pattern_long <- pattern_long %>%
    dplyr::left_join(pattern_df %>% dplyr::select(pattern_id, x_pos), by = "pattern_id")

  # Set up intersection bar colors based on highlighting
  intersection_bar_colors <- rep(main_bar_color, nrow(pattern_df))
  names(intersection_bar_colors) <- pattern_df$pattern_id

  # Apply highlights if specified
  highlight_indices <- NULL
  if (!is.null(highlight_intersections)) {
    if (is.list(highlight_intersections) || is.vector(highlight_intersections)) {
      valid_highlights <- highlight_intersections[highlight_intersections %in% pattern_df$pattern_id]
      intersection_bar_colors[as.character(valid_highlights)] <- highlight_color
      highlight_indices <- valid_highlights
    }
  }

  # Create a minimal theme for the plots
  custom_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )

  # Define limits for alignment
  x_limits <- c(0.5, nrow(pattern_df) + 0.5)
  y_limits <- c(0.5, nrow(set_size_df) + 0.5)

  # Apply highlight colors to dots if highlight_indices is provided
  if (!is.null(highlight_indices)) {
    pattern_long$is_highlighted <- pattern_long$pattern_id %in% highlight_indices
  } else {
    pattern_long$is_highlighted <- FALSE
  }

  # Create a dot color vector based on highlight status
  pattern_long$dot_color <- ifelse(pattern_long$is_in_set == 1,
                                   ifelse(pattern_long$is_highlighted, highlight_color, intersection_color),
                                   "white")

  # Create the matrix visualization
  matrix_plot <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = pattern_long,
      ggplot2::aes(
        x = x_pos,
        y = y_pos,
        size = ifelse(is_in_set == 1, point_size, point_size/2)
      ),
      shape = 21,
      fill = pattern_long$dot_color,
      color = "black",
      stroke = 0.3
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::labs(x = "", y = "") +
    custom_theme +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = text_size, angle = text_angle),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 0, r = 0.5, b = 0.5, l = 0, unit = "cm")
    ) +
    ggplot2::scale_x_continuous(limits = x_limits, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      breaks = set_size_df$y_pos,
      labels = set_size_df$set_name,
      limits = y_limits,
      expand = c(0, 0)
    )

  # Add connecting lines for intersections
  # Group data by pattern_id to easily extract sets for each intersection
  grouped_data <- pattern_long %>%
    dplyr::filter(is_in_set == 1) %>%
    dplyr::group_by(pattern_id, x_pos, is_highlighted) %>%
    dplyr::summarise(
      sets = list(y_pos),
      n_sets = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_sets > 1)  # Only keep intersections with multiple sets

  # Add lines for each intersection with multiple sets
  for (i in 1:nrow(grouped_data)) {
    # Extract y positions
    y_positions <- unlist(grouped_data$sets[i])

    line_data <- data.frame(
      x_pos = grouped_data$x_pos[i],
      y_pos = y_positions,
      is_highlighted = grouped_data$is_highlighted[i]
    )

    # Line color based on highlight status
    line_color <- if (line_data$is_highlighted[1]) highlight_color else intersection_color

    # Add the line
    matrix_plot <- matrix_plot +
      ggplot2::geom_line(
        data = line_data,
        ggplot2::aes(x = x_pos, y = y_pos),
        size = line_size,
        color = line_color
      )
  }

  # Create the intersection size bar plot (top)
  intersection_plot <- ggplot2::ggplot() +
    ggplot2::geom_col(
      data = pattern_df,
      ggplot2::aes(
        x = x_pos,
        y = size,
        fill = factor(pattern_id)
      ),
      width = bar_width
    ) +
    ggplot2::scale_fill_manual(values = stats::setNames(intersection_bar_colors, 1:length(intersection_bar_colors))) +
    ggplot2::labs(x = "", y = intersection_title) +
    custom_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 0.5, r = 0.5, b = 0, l = 0.5, unit = "cm")
    ) +
    ggplot2::scale_x_continuous(limits = x_limits, expand = c(0, 0))

  # Add size labels if requested
  if (show_numbers) {
    intersection_plot <- intersection_plot +
      ggplot2::geom_text(
        data = pattern_df,
        ggplot2::aes(
          x = x_pos,
          y = size,
          label = size
        ),
        vjust = -0.5,
        size = 3
      )
  }

  # Create the set size bar plot (left)
  set_size_plot <- ggplot2::ggplot() +
    ggplot2::geom_col(
      data = set_size_df,
      ggplot2::aes(
        y = y_pos,
        x = size,
        fill = set_name
      ),
      width = bar_width
    ) +
    ggplot2::labs(x = set_size_title, y = "") +
    ggplot2::scale_fill_manual(values = sets_bar_colors) +
    custom_theme +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 0, r = 0, b = 0.5, l = 0.5, unit = "cm")
    ) +
    ggplot2::scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
    ggplot2::scale_x_continuous(position = "top")

  # Add size labels to set size plot
  set_size_plot <- set_size_plot +
    ggplot2::geom_text(
      data = set_size_df,
      ggplot2::aes(x = size/2, y = y_pos, label = size),
      size = 3
    )

  # Create an empty plot for the top-left corner
  empty_plot <- patchwork::plot_spacer()

  # Use patchwork for layout management
  # First create the top row
  top_row <- empty_plot + intersection_plot +
    patchwork::plot_layout(widths = c(0.3, 0.7))

  # Then create the bottom row
  bottom_row <- set_size_plot + matrix_plot +
    patchwork::plot_layout(widths = c(0.3, 0.7))

  # Combine the rows
  combined_plot <- top_row / bottom_row +
    patchwork::plot_layout(heights = c(0.4, 0.6))

  # Add title if specified
  if (!is.null(title)) {
    combined_plot <- combined_plot + patchwork::plot_annotation(title = title)
  }

  return(combined_plot)
}

#' Create an UpSet Plot for DEG Analysis
#'
#' Creates an UpSet plot to visualize shared and unique differentially expressed genes
#' across different cell types or comparisons.
#'
#' @param deg_results Output from findDifferentialGenes or compareDEGs function
#' @param by_cell_type If TRUE, shows overlap across cell types; if FALSE, across comparisons
#' @param direction DEG direction to include: "up", "down", or "both"
#' @param min_size Minimum intersection size to include
#' @param cell_types Vector of cell types to include (default: all)
#' @param comparisons Vector of comparisons to include (default: all)
#' @param use_gene_id If TRUE, uses gene IDs instead of symbols when available
#' @param point_size Size of points in the matrix
#' @param main_bar_color Color for the intersection size bars
#' @param title Plot title
#' @param show_numbers Whether to show numbers on bars
#' @param ... Additional parameters passed to upsetPlot function
#'
#' @return An UpSet plot (patchwork object)
#' @export
visualizeDEGUpset <- function(deg_results,
                              by_cell_type = TRUE,
                              direction = "both",
                              min_size = 5,
                              cell_types = NULL,
                              comparisons = NULL,
                              use_gene_id = FALSE,
                              point_size = 3,
                              main_bar_color = "steelblue",
                              title = NULL,
                              show_numbers = TRUE,
                              ...) {

  # First, prepare the data based on input type
  if (inherits(deg_results, "list")) {
    # Extract DEG results from object
    if ("all_degs" %in% names(deg_results)) {
      all_degs <- deg_results$all_degs
    } else if ("degs" %in% names(deg_results)) {
      # Combine individual result tables
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

  # Ensure required columns exist
  req_cols <- c("gene", "significant")
  if (!all(req_cols %in% colnames(all_degs))) {
    stop(paste("Missing required columns:", paste(req_cols[!req_cols %in% colnames(all_degs)], collapse = ", ")))
  }

  # Filter for significant genes only
  sig_degs <- all_degs[all_degs$significant, ]

  # Filter by direction if specified
  if (direction == "up") {
    if ("direction" %in% colnames(sig_degs)) {
      sig_degs <- sig_degs[sig_degs$direction == "up", ]
    } else if ("avg_log2FC" %in% colnames(sig_degs)) {
      sig_degs <- sig_degs[sig_degs$avg_log2FC > 0, ]
    } else {
      warning("Cannot filter by direction - no direction or avg_log2FC column found")
    }
  } else if (direction == "down") {
    if ("direction" %in% colnames(sig_degs)) {
      sig_degs <- sig_degs[sig_degs$direction == "down", ]
    } else if ("avg_log2FC" %in% colnames(sig_degs)) {
      sig_degs <- sig_degs[sig_degs$avg_log2FC < 0, ]
    } else {
      warning("Cannot filter by direction - no direction or avg_log2FC column found")
    }
  }

  # Check if we still have data after filtering
  if (nrow(sig_degs) == 0) {
    stop("No significant DEGs found after filtering")
  }

  # Filter by cell types if specified
  if (!is.null(cell_types)) {
    if (!all(cell_types %in% unique(sig_degs$cell_type))) {
      warning("Some cell types not found in the data")
    }
    sig_degs <- sig_degs[sig_degs$cell_type %in% cell_types, ]
  }

  # Filter by comparisons if specified
  if (!is.null(comparisons)) {
    if (!all(comparisons %in% unique(sig_degs$comparison))) {
      warning("Some comparisons not found in the data")
    }
    sig_degs <- sig_degs[sig_degs$comparison %in% comparisons, ]
  }

  # Use gene ID instead of symbol if requested and available
  gene_col <- "gene"
  if (use_gene_id && "id" %in% colnames(sig_degs)) {
    gene_col <- "id"
  }

  # Create gene lists based on grouping
  gene_lists <- list()

  if (by_cell_type) {
    # Create lists of genes by cell type
    for (ct in unique(sig_degs$cell_type)) {
      ct_genes <- unique(sig_degs[sig_degs$cell_type == ct, gene_col])
      if (length(ct_genes) > 0) {
        gene_lists[[ct]] <- ct_genes
      }
    }

    # Set title if not provided
    if (is.null(title)) {
      title_direction <- ifelse(direction == "both", "DEGs",
                                ifelse(direction == "up", "Up-regulated Genes", "Down-regulated Genes"))
      title <- paste(title_direction, "Across Cell Types")
    }
  } else {
    # Create lists of genes by comparison
    for (comp in unique(sig_degs$comparison)) {
      comp_genes <- unique(sig_degs[sig_degs$comparison == comp, gene_col])
      if (length(comp_genes) > 0) {
        gene_lists[[comp]] <- comp_genes
      }
    }

    # Set title if not provided
    if (is.null(title)) {
      title_direction <- ifelse(direction == "both", "DEGs",
                                ifelse(direction == "up", "Up-regulated Genes", "Down-regulated Genes"))
      title <- paste(title_direction, "Across Comparisons")
    }
  }

  # Create the upset plot
  upset_plot <- upsetPlot(
    data_list = gene_lists,
    min_intersection_size = min_size,
    main_bar_color = main_bar_color,
    point_size = point_size,
    title = title,
    show_numbers = show_numbers,
    ...
  )

  return(upset_plot)
}
