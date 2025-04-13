#' Create an Enhanced Gene Regulation Heatmap
#'
#' Creates a heatmap visualizing the up and down regulated genes with
#' row clustering, regulation indicators, and labels positioned properly above the heatmap.
#'
#' @param data A data frame containing DEG results
#' @param comparison_col Column name for the comparison categories
#' @param cell_type_col Column name for the cell types
#' @param up_col Column name for up-regulated gene counts
#' @param down_col Column name for down-regulated gene counts
#' @param border_color Color for tile borders (default = "white")
#' @param border_size Size of tile borders (default = 0.5)
#' @param show_counts Logical. Whether to show count numbers inside tiles (default = TRUE)
#' @param count_size Size of count numbers (default = 3)
#' @param count_color Color of count numbers (default = "black")
#' @param low_color Color for low values in gradient (default = "deepskyblue")
#' @param mid_color Color for mid values in gradient (default = "white")
#' @param high_color Color for high values in gradient (default = "darkorange")
#' @param cluster_rows Logical. Whether to cluster rows (default = TRUE)
#' @param cluster_cols Logical. Whether to cluster columns (default = FALSE)
#' @param cluster_method Clustering method (default = "complete")
#' @param distance_method Distance method for clustering (default = "euclidean")
#' @param group_regulation Logical. If TRUE, groups all up and all down together (default = FALSE)
#' @param show_reg_indicators Logical. If TRUE, adds color indicators for up/down regulation (default = TRUE)
#' @param up_indicator_color Color for up-regulation indicator (default = "darkorange")
#' @param down_indicator_color Color for down-regulation indicator (default = "deepskyblue")
#' @param indicator_height Height of the color indicator (default = 0.2)
#' @param top_labels Logical. If TRUE, shows Up/Down labels at top of heatmap (default = TRUE)
#' @param top_label_size Size of the top labels (default = 3.5)
#' @param label_distance Distance of labels from the top of the heatmap (default = 1.5)
#' @param title Plot title (default = "Up and Down Regulated Genes by Cell Type and Comparison")
#' @param title_size Size of the title (default = 12)
#' @param x_lab Label for x-axis (default = "Comparison")
#' @param y_lab Label for y-axis (default = "Cell Type")
#' @param legend_title Title for the legend (default = "Gene Count")
#' @param x_angle Angle for x-axis text (default = 45)
#' @param top_margin Extra space at the top of the plot (default = 40)
#'
#' @return A ggplot object
#'
gene_regulation_heatmap_enhanced <- function(data,
                                            comparison_col = "comparison_pattern",
                                            cell_type_col = "cell_type",
                                            up_col = "up_genes",
                                            down_col = "down_genes",
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
                                            title = "DEGs Across Across Cell Type and Comparisons",
                                            title_size = 12,
                                            x_lab = "Comparison",
                                            y_lab = "Cell Type",
                                            legend_title = "Gene Count",
                                            x_angle = 45,
                                            top_margin = 40) {
  
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(scales)
  require(stats)
  
  # Create a copy of the data to avoid modifying the original
  plot_data <- data
  
  # Validate column names
  col_names <- colnames(plot_data)
  if (!comparison_col %in% col_names) stop(paste("Column", comparison_col, "not found in data"))
  if (!cell_type_col %in% col_names) stop(paste("Column", cell_type_col, "not found in data"))
  if (!up_col %in% col_names) stop(paste("Column", up_col, "not found in data"))
  if (!down_col %in% col_names) stop(paste("Column", down_col, "not found in data"))
  
  # Create long format data
  side_data <- plot_data %>%
    pivot_longer(
      cols = c(up_col, down_col),
      names_to = "regulation",
      values_to = "count"
    ) %>%
    mutate(regulation = ifelse(regulation == up_col, "Up", "Down"))
  
  # Reshape data to matrix format for clustering
  wide_data <- side_data %>%
    unite("comparison_regulation", c(!!sym(comparison_col), regulation), sep = "_") %>%
    pivot_wider(
      id_cols = !!sym(cell_type_col),
      names_from = comparison_regulation,
      values_from = count,
      values_fill = 0
    )
  
  # Extract row names (cell types)
  cell_types <- wide_data[[cell_type_col]]
  
  # Convert to matrix for clustering, excluding cell_type column
  mat_data <- as.matrix(wide_data[, -1])
  rownames(mat_data) <- cell_types
  
  # Cluster rows if requested
  if (cluster_rows) {
    row_dist <- dist(mat_data, method = distance_method)
    row_clust <- hclust(row_dist, method = cluster_method)
    cell_type_order <- rownames(mat_data)[row_clust$order]
    side_data[[cell_type_col]] <- factor(side_data[[cell_type_col]], levels = cell_type_order)
  } else {
    side_data[[cell_type_col]] <- factor(side_data[[cell_type_col]])
  }
  
  # Handle comparison ordering
  comparisons <- unique(side_data[[comparison_col]])
  
  if (cluster_cols) {
    # Create matrix with comparison averages for clustering
    comparison_names <- comparisons
    comp_matrix <- matrix(0, nrow = length(cell_types), ncol = length(comparison_names))
    rownames(comp_matrix) <- cell_types
    colnames(comp_matrix) <- comparison_names
    
    for (i in 1:length(comparison_names)) {
      comp <- comparison_names[i]
      up_cols <- grepl(paste0(comp, "_Up"), colnames(mat_data))
      down_cols <- grepl(paste0(comp, "_Down"), colnames(mat_data))
      comp_matrix[, i] <- rowMeans(cbind(mat_data[, up_cols, drop = FALSE], 
                                         mat_data[, down_cols, drop = FALSE]))
    }
    
    col_dist <- dist(t(comp_matrix), method = distance_method)
    col_clust <- hclust(col_dist, method = cluster_method)
    comparison_order <- colnames(comp_matrix)[col_clust$order]
    side_data[[comparison_col]] <- factor(side_data[[comparison_col]], levels = comparison_order)
  } else {
    # If not clustering, still make comparison a factor to preserve order
    side_data[[comparison_col]] <- factor(side_data[[comparison_col]])
  }
  
  # Handle regulation grouping
  if (group_regulation) {
    # Create a combined x-axis position that groups all up together and all down together
    side_data <- side_data %>%
      mutate(
        x_pos = ifelse(regulation == "Up",
                      as.numeric(factor(!!sym(comparison_col))),
                      as.numeric(factor(!!sym(comparison_col))) + length(levels(factor(!!sym(comparison_col)))) + 1)
      )
    
    # Create positions for comparison labels
    up_label_positions <- 1:length(comparisons)
    down_label_positions <- (length(comparisons) + 2):(2*length(comparisons) + 1)
    all_label_positions <- c(up_label_positions, down_label_positions)
    
    # Create labels with line breaks for readability
    all_labels <- rep(levels(factor(side_data[[comparison_col]])), 2)
    
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
    p <- ggplot() +
      # Add background rectangles for up/down groups if indicators are requested
      {if(show_reg_indicators) geom_rect(
        data = indicator_df,
        aes(xmin = min(up_label_positions) - 0.5, 
            xmax = max(up_label_positions) + 0.5,
            ymin = -Inf, 
            ymax = Inf),
        fill = up_indicator_color,
        alpha = 0.1
      )} +
      {if(show_reg_indicators) geom_rect(
        data = indicator_df,
        aes(xmin = min(down_label_positions) - 0.5, 
            xmax = max(down_label_positions) + 0.5,
            ymin = -Inf, 
            ymax = Inf),
        fill = down_indicator_color,
        alpha = 0.1
      )} +
      # Add the main tiles
      geom_tile(
        data = side_data,
        aes(x = x_pos, y = !!sym(cell_type_col), fill = count),
        color = border_color, 
        size = border_size, 
        width = 0.9, 
        height = 0.9
      ) +
      # Add count numbers if requested
      {if(show_counts) geom_text(
        data = side_data,
        aes(x = x_pos, y = !!sym(cell_type_col), label = count),
        size = count_size,
        color = count_color
      )} +
      # Custom color scale
      scale_fill_gradientn(
        colors = colorRampPalette(c(low_color, mid_color, high_color))(100),
        name = legend_title
      ) +
      # Set x-axis with custom breaks and labels
      scale_x_continuous(
        breaks = x_breaks,
        labels = x_labels
      )
    
    # Add group labels at the top if requested
    if (top_labels) {
      max_y <- max(as.numeric(factor(side_data[[cell_type_col]])))
      p <- p + 
        annotate("text", 
                x = mean(up_label_positions), 
                y = max_y + label_distance, 
                label = "Up-regulated", 
                size = top_label_size, 
                fontface = "bold",
                color = up_indicator_color) +
        annotate("text", 
                x = mean(down_label_positions), 
                y = max_y + label_distance, 
                label = "Down-regulated", 
                size = top_label_size, 
                fontface = "bold",
                color = down_indicator_color)
    }
      
  } else {
    # Side-by-side approach (original)
    side_data <- side_data %>%
      mutate(
        comp_id = as.numeric(factor(!!sym(comparison_col), levels = levels(factor(!!sym(comparison_col))))),
        x_pos = comp_id * 2 + ifelse(regulation == "Up", -0.5, 0.5)
      )
    
    # Create secondary x-axis labels for Up/Down
    comp_ids <- unique(side_data$comp_id)
    up_positions <- comp_ids * 2 - 0.5
    down_positions <- comp_ids * 2 + 0.5
    
    # Calculate the midpoints for the up and down sections
    up_mid <- mean(range(up_positions))
    down_mid <- mean(range(down_positions))
    
    # Create the plot
    p <- ggplot()
    
    # Add regulation indicators at the top if requested
    if (show_reg_indicators) {
      # Define the maximum y position for consistent plotting
      max_y_level <- length(levels(factor(side_data[[cell_type_col]])))
      
      # Add color indicators at the top of the plot for Up regulation
      for (pos in up_positions) {
        p <- p + 
          annotate("rect", 
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
          annotate("rect", 
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
          annotate("text", 
                 x = up_mid, 
                 y = max_y_level + label_distance, 
                 label = "Up-regulated", 
                 size = top_label_size,
                 fontface = "bold",
                 color = up_indicator_color) +
          annotate("text", 
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
      geom_point(
        data = dummy_data,
        aes(x = x, y = y, color = regulation),
        shape = 15,   # Square shape for better legend appearance
        size = 5,     # Larger size for visibility in legend
        alpha = 0     # Make invisible in the plot
      ) +
      # Define regulation colors
      scale_color_manual(
        values = reg_colors,
        name = "Regulation"
      )
    
    # Add the main tiles
    p <- p + 
      geom_tile(
        data = side_data,
        aes(x = x_pos, y = !!sym(cell_type_col), fill = count),
        color = border_color, 
        size = border_size, 
        width = 0.9, 
        height = 0.9
      )
    
    # Add count numbers if requested
    if (show_counts) {
      p <- p + geom_text(
        data = side_data,
        aes(x = x_pos, y = !!sym(cell_type_col), label = count),
        size = count_size,
        color = count_color
      )
    }
    
    # Custom color scale
    p <- p + scale_fill_gradientn(
      colors = colorRampPalette(c(low_color, mid_color, high_color))(100),
      name = legend_title
    )
    
    # Set x-axis with labels for comparisons (no secondary axis)
    p <- p + scale_x_continuous(
      breaks = comp_ids * 2,
      labels = levels(factor(side_data[[comparison_col]]))
    )
  }
  
  # Define the maximum y position for expansion, based on label distance
  max_y_expansion <- if(top_labels) label_distance + 0.5 else indicator_height + 0.3
  
  # Complete the plot styling
  p <- p + 
    labs(
      title = title,
      x = x_lab,
      y = y_lab
    ) +
    # Expand y-axis to make room for the top labels
    scale_y_discrete(
      expand = expansion(add = c(0, max_y_expansion))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = title_size, face = "bold"),
      axis.text.x = element_text(angle = x_angle, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = top_margin, r = 20, b = 20, l = 20, unit = "pt"),
      legend.position = "right"
    )
  
  return(p)
}

# Example usage with adjusted label positioning:
# gene_regulation_heatmap_enhanced(
#   deg_results, 
#   show_reg_indicators = TRUE,
#   top_labels = TRUE,
#   label_distance = 1.5,      # Increased distance from heatmap
#   top_label_size = 3.5,      # Slightly larger labels
#   top_margin = 40           # Increased top margin
# )