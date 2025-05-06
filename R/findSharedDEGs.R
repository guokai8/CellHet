# findSharedDEGs.R

#' Analyze Shared and Unique DEGs Across Comparisons
#'
#' @param deg_results Output from compareDEGs function
#' @param cell_types Vector of cell types to include (default: all)
#' @param min_deg_count Minimum number of DEGs required to include a comparison
#' @param plot_type Type of plot ("upset", "venn", "heatmap", "network")
#' @param direction Filter for gene direction ("up", "down", "both")
#'
#' @return A list with overlap analysis and corresponding visualization
#' @export
findSharedDEGs <- function(deg_results,
                           cell_types = NULL,
                           min_deg_count = 10,
                           plot_type = "upset",
                           direction = "both") {

  # Extract DEG results
  deg_list <- deg_results$degs

  if (length(deg_list) == 0) {
    stop("No DEG results found in the input data")
  }

  # Filter by cell type if specified
  if (!is.null(cell_types)) {
    # Extract cell types from result keys
    result_cell_types <- sapply(names(deg_list), function(x) strsplit(x, "__")[[1]][1])
    deg_list <- deg_list[result_cell_types %in% cell_types]

    if (length(deg_list) == 0) {
      stop("No results found for specified cell types")
    }
  }

  # Extract significant genes for each comparison
  gene_sets <- list()
  for (comparison_name in names(deg_list)) {
    deg_data <- deg_list[[comparison_name]]

    # Filter by direction and significance
    if (direction == "up") {
      genes <- deg_data$gene[deg_data$significant & deg_data$avg_log2FC > 0]
    } else if (direction == "down") {
      genes <- deg_data$gene[deg_data$significant & deg_data$avg_log2FC < 0]
    } else {
      genes <- deg_data$gene[deg_data$significant]
    }

    # Only include if enough DEGs
    if (length(genes) >= min_deg_count) {
      gene_sets[[comparison_name]] <- genes
    }
  }

  if (length(gene_sets) == 0) {
    stop("No comparisons with sufficient DEGs after filtering")
  }

  # Find all unique genes
  all_genes <- unique(unlist(gene_sets))

  # Create presence/absence matrix
  presence_matrix <- matrix(0, nrow = length(all_genes), ncol = length(gene_sets))
  rownames(presence_matrix) <- all_genes
  colnames(presence_matrix) <- names(gene_sets)

  for (i in seq_along(gene_sets)) {
    presence_matrix[gene_sets[[i]], i] <- 1
  }

  # Calculate gene set overlaps
  overlaps <- list()
  all_pairs <- combn(names(gene_sets), 2, simplify = FALSE)

  for (pair in all_pairs) {
    set1 <- gene_sets[[pair[1]]]
    set2 <- gene_sets[[pair[2]]]
    common_genes <- intersect(set1, set2)

    overlaps[[paste(pair, collapse = "_vs_")]] <- list(
      set1_name = pair[1],
      set2_name = pair[2],
      common_genes = common_genes,
      set1_only = setdiff(set1, set2),
      set2_only = setdiff(set2, set1),
      jaccard_index = length(common_genes) / length(union(set1, set2))
    )
  }

  # Create visualization based on plot_type
  if (plot_type == "upset") {
    # Use our custom upset plot implementation
    plot <- visualizeDEGUpset(
      deg_results = deg_results,
      by_cell_type = TRUE,
      direction = direction,
      min_size = min_deg_count,
      cell_types = cell_types
    )
  } else if (plot_type == "venn") {
    if (length(gene_sets) > 5) {
      warning("Venn diagrams not recommended for more than 5 sets. Consider using 'upset' instead.")
    }

    if (!requireNamespace("VennDiagram", quietly = TRUE)) {
      stop("Package 'VennDiagram' needed for venn plots. Please install it.")
    }

    # Create temporary file for venn diagram
    temp_file <- tempfile(fileext = ".png")

    # Generate venn diagram
    venn_plot <- VennDiagram::venn.diagram(
      x = gene_sets,
      filename = temp_file,
      fill = rainbow(length(gene_sets)),
      alpha = 0.5,
      cex = 1,
      cat.cex = 1,
      cat.fontface = "bold"
    )

    # Read the image back
    if (!requireNamespace("png", quietly = TRUE)) {
      warning("Package 'png' needed to display venn diagram. Returning NULL plot.")
      plot <- NULL
    } else {
      plot <- png::readPNG(temp_file)
      # Delete temporary file
      unlink(temp_file)
    }

  } else if (plot_type == "heatmap") {
    # Create heatmap of shared genes
    gene_counts <- rowSums(presence_matrix)
    shared_genes <- names(gene_counts[gene_counts > 1])

    if (length(shared_genes) == 0) {
      warning("No shared genes found. Consider using 'min_deg_count' parameter.")
      plot <- NULL
    } else {
      # Subset matrix for visualization
      plot_matrix <- presence_matrix[shared_genes, ]

      # Create heatmap
      if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        # Use ComplexHeatmap
        plot <- ComplexHeatmap::Heatmap(
          matrix = plot_matrix,
          name = "Present",
          row_title = "Genes",
          column_title = "Comparisons",
          row_names_gp = grid::gpar(fontsize = 10),
          column_names_gp = grid::gpar(fontsize = 10),
          column_names_rot = 45
        )
      } else {
        # Fallback to geom_tile
        heatmap_data <- reshape2::melt(plot_matrix,
                                       varnames = c("gene", "comparison"),
                                       value.name = "present")

        plot <- ggplot(heatmap_data, aes(x = comparison, y = gene, fill = factor(present))) +
          geom_tile() +
          scale_fill_manual(values = c("0" = "white", "1" = "steelblue")) +
          theme_minimal() +
          labs(
            title = "Shared DEGs across comparisons",
            x = "Comparison",
            y = "Gene",
            fill = "Present"
          ) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
    }

  } else if (plot_type == "network") {
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("Package 'igraph' needed for network plots. Please install it.")
    }

    # Create network from overlaps
    edges <- do.call(rbind, lapply(names(overlaps), function(name) {
      overlap <- overlaps[[name]]
      data.frame(
        from = overlap$set1_name,
        to = overlap$set2_name,
        weight = length(overlap$common_genes),
        jaccard = overlap$jaccard_index
      )
    }))

    # Create graph
    g <- igraph::graph_from_data_frame(edges, directed = FALSE,
                                       vertices = data.frame(name = names(gene_sets)))

    # Add vertex attributes
    igraph::V(g)$size <- sapply(gene_sets, length)

    # Add edge attributes
    igraph::E(g)$width <- sqrt(edges$weight)

    # Create network plot
    if (requireNamespace("ggraph", quietly = TRUE)) {
      # Use ggraph for better visualization
      plot <- ggraph::ggraph(g, layout = "fr") +
        ggraph::geom_edge_link(aes(width = width, alpha = jaccard),
                               edge_colour = "skyblue") +
        ggraph::geom_node_point(aes(size = size), color = "steelblue") +
        ggraph::geom_node_text(aes(label = name), repel = TRUE) +
        ggraph::scale_edge_width(range = c(0.5, 3)) +
        ggraph::scale_edge_alpha(range = c(0.2, 1)) +
        ggraph::theme_graph() +
        labs(title = "DEG Sharing Network")
    } else {
      # Fallback to basic igraph plot
      plot <- igraph::plot.igraph(
        g,
        vertex.size = sqrt(igraph::V(g)$size) * 0.1,
        vertex.label = igraph::V(g)$name,
        edge.width = igraph::E(g)$width * 0.5,
        layout = igraph::layout_with_fr(g)
      )
    }
  } else {
    stop(paste0("Unsupported plot type: ", plot_type))
  }

  # Return results
  return(list(
    gene_sets = gene_sets,
    presence_matrix = presence_matrix,
    overlaps = overlaps,
    plot = plot
  ))
}
