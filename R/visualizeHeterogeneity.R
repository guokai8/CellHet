# visualizeHeterogeneity.R

#' Create Heterogeneity Map Across Cell Types
#'
#' @param deg_results Output from compareDEGs function
#' @param pathway_results Optional output from runPathwayAnalysis
#' @param layout_type Visualization type ("heatmap", "network", "umap", "circular")
#' @param level Analysis level ("gene", "pathway", "both")
#' @param similarity_metric Metric to calculate similarity between cell types
#'
#' @return A plot object showing heterogeneity across cell types
#' @export
visualizeHeterogeneity <- function(deg_results,
                                   pathway_results = NULL,
                                   layout_type = "heatmap",
                                   level = "gene",
                                   similarity_metric = "jaccard") {

  # Check input
  if (level == "pathway" && is.null(pathway_results)) {
    stop("Pathway results required for pathway-level heterogeneity analysis")
  }

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

  # Create similarity matrix between cell types
  similarity_matrix <- matrix(0,
                              nrow = length(unique_cell_types),
                              ncol = length(unique_cell_types))
  rownames(similarity_matrix) <- unique_cell_types
  colnames(similarity_matrix) <- unique_cell_types

  # Fill similarity matrix based on specified level and metric
  for (i in 1:length(unique_cell_types)) {
    for (j in i:length(unique_cell_types)) {
      ct1 <- unique_cell_types[i]
      ct2 <- unique_cell_types[j]

      if (i == j) {
        # Same cell type, perfect similarity
        similarity_matrix[i, j] <- 1
      } else {
        # Calculate similarity between cell types
        if (level == "gene" || level == "both") {
          # For each comparison, get DEGs for both cell types
          gene_similarity <- calculateCellTypeSimilarity(
            deg_list, ct1, ct2, unique_comparisons, "gene", similarity_metric
          )

          if (level == "gene") {
            similarity_matrix[i, j] <- similarity_matrix[j, i] <- gene_similarity
          }
        }

        if (level == "pathway" || level == "both") {
          # For each comparison, get pathways for both cell types
          pathway_similarity <- calculateCellTypeSimilarity(
            pathway_results$pathway_results, ct1, ct2, unique_comparisons, "pathway", similarity_metric
          )

          if (level == "pathway") {
            similarity_matrix[i, j] <- similarity_matrix[j, i] <- pathway_similarity
          } else if (level == "both") {
            # Average gene and pathway similarity
            similarity_matrix[i, j] <- similarity_matrix[j, i] <- (gene_similarity + pathway_similarity) / 2
          }
        }
      }
    }
  }

  # Create heterogeneity score (1 - similarity)
  heterogeneity_matrix <- 1 - similarity_matrix

  # Create visualization based on layout_type
  if (layout_type == "heatmap") {
    # Create heatmap of heterogeneity
    if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      # Use ComplexHeatmap
      plot <- ComplexHeatmap::Heatmap(
        matrix = heterogeneity_matrix,
        name = "Heterogeneity",
        row_title = "Cell Type",
        column_title = "Cell Type",
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 10),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid::grid.text(sprintf("%.2f", heterogeneity_matrix[i, j]), x, y,
                          gp = grid::gpar(fontsize = 10))
        }
      )
    } else {
      # Fallback to geom_tile
      heatmap_data <- reshape2::melt(heterogeneity_matrix,
                                     varnames = c("cell_type1", "cell_type2"),
                                     value.name = "heterogeneity")

      plot <- ggplot(heatmap_data,
                     aes(x = cell_type1, y = cell_type2, fill = heterogeneity)) +
        geom_tile() +
        geom_text(aes(label = sprintf("%.2f", heterogeneity))) +
        scale_fill_viridis_c(option = "magma") +
        theme_minimal() +
        labs(
          title = "Cell Type Heterogeneity Map",
          x = "Cell Type",
          y = "Cell Type",
          fill = "Heterogeneity"
        )
    }

  } else if (layout_type == "network") {
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("Package 'igraph' needed for network plots. Please install it.")
    }

    # Create network from heterogeneity matrix
    # Lower heterogeneity = stronger edge
    edge_strength <- 1 - heterogeneity_matrix
    edge_strength[edge_strength < 0.2] <- 0  # Threshold to remove weak edges

    # Create graph
    g <- igraph::graph_from_adjacency_matrix(
      edge_strength,
      mode = "undirected",
      weighted = TRUE
    )

    # Community detection for coloring
    comm <- igraph::cluster_louvain(g)
    igraph::V(g)$community <- comm$membership

    # Calculate node sizes based on DEG counts
    deg_counts <- sapply(unique_cell_types, function(ct) {
      # Get indices for this cell type
      idx <- which(cell_types == ct)
      # Sum DEGs across all comparisons
      sum(sapply(deg_list[result_keys[idx]], function(deg) sum(deg$significant)))
    })

    igraph::V(g)$size <- log10(deg_counts + 1) * 5

    # Create network plot
    if (requireNamespace("ggraph", quietly = TRUE)) {
      # Use ggraph for better visualization
      plot <- ggraph::ggraph(g, layout = "fr") +
        ggraph::geom_edge_link(aes(width = weight, alpha = weight),
                               edge_colour = "grey50") +
        ggraph::geom_node_point(aes(size = size, color = factor(community))) +
        ggraph::geom_node_text(aes(label = name), repel = TRUE) +
        ggraph::scale_edge_width(range = c(0.1, 2)) +
        ggraph::scale_edge_alpha(range = c(0.1, 1)) +
        ggraph::theme_graph() +
        labs(title = "Cell Type Heterogeneity Network",
             subtitle = paste0("Based on ", level, " level ", similarity_metric, " similarity"))
    } else {
      # Fallback to basic igraph plot
      plot <- igraph::plot.igraph(
        g,
        vertex.size = igraph::V(g)$size,
        vertex.label = igraph::V(g)$name,
        vertex.color = igraph::V(g)$community,
        edge.width = igraph::E(g)$weight * 3,
        layout = igraph::layout_with_fr(g)
      )
    }

  } else if (layout_type == "umap") {
    if (!requireNamespace("umap", quietly = TRUE)) {
      stop("Package 'umap' needed for UMAP plots. Please install it.")
    }

    # Run UMAP on the heterogeneity matrix
    umap_result <- umap::umap(heterogeneity_matrix)

    # Create data frame for plotting
    umap_df <- data.frame(
      cell_type = unique_cell_types,
      x = umap_result$layout[, 1],
      y = umap_result$layout[, 2]
    )

    # Calculate DEG counts for point sizes
    umap_df$deg_count <- sapply(unique_cell_types, function(ct) {
      # Get indices for this cell type
      idx <- which(cell_types == ct)
      # Sum DEGs across all comparisons
      sum(sapply(deg_list[result_keys[idx]], function(deg) sum(deg$significant)))
    })

    # Create UMAP plot
    plot <- ggplot(umap_df, aes(x = x, y = y)) +
      geom_point(aes(size = log10(deg_count + 1), color = cell_type), alpha = 0.7) +
      geom_text_repel(aes(label = cell_type), size = 4, box.padding = 0.5) +
      theme_minimal() +
      labs(
        title = "Cell Type Heterogeneity UMAP",
        subtitle = paste0("Based on ", level, " level ", similarity_metric, " similarity"),
        x = "UMAP1",
        y = "UMAP2",
        size = "log10(DEGs)",
        color = "Cell Type"
      ) +
      theme(legend.position = "right")

  } else if (layout_type == "circular") {
    if (!requireNamespace("circlize", quietly = TRUE)) {
      stop("Package 'circlize' needed for circular plots. Please install it.")
    }

    # Create data for chord diagram
    # Use heterogeneity threshold to determine connections
    threshold <- 0.5
    chord_matrix <- heterogeneity_matrix
    chord_matrix[chord_matrix > threshold] <- 0

    # Normalize for better visualization
    chord_matrix <- chord_matrix / max(chord_matrix)

    # Plot chord diagram
    circlize::circos.clear()
    circlize::circos.par(gap.degree = 5)

    # Generate color for each cell type
    cell_type_colors <- rainbow(length(unique_cell_types))
    names(cell_type_colors) <- unique_cell_types

    # Create plot
    plot <- circlize::chordDiagram(
      chord_matrix,
      grid.col = cell_type_colors,
      transparency = 0.2,
      directional = FALSE
    )

    # Add title
    title(main = "Cell Type Heterogeneity Chord Diagram",
          sub = paste0("Based on ", level, " level ", similarity_metric, " similarity"))

  } else {
    stop(paste0("Unsupported layout type: ", layout_type))
  }

  # Return results
  return(list(
    plot = plot,
    heterogeneity_matrix = heterogeneity_matrix,
    similarity_matrix = similarity_matrix
  ))
}

#' Calculate Similarity Between Two Cell Types
#'
#' @param result_list List of results (DEGs or pathways)
#' @param cell_type1 First cell type
#' @param cell_type2 Second cell type
#' @param comparisons Vector of comparison names
#' @param level Analysis level ("gene" or "pathway")
#' @param method Similarity metric
#'
#' @return Numeric similarity score
#' @keywords internal
calculateCellTypeSimilarity <- function(result_list,
                                        cell_type1,
                                        cell_type2,
                                        comparisons,
                                        level = "gene",
                                        method = "jaccard") {

  # Initialize similarities for each comparison
  comparison_similarities <- numeric(length(comparisons))

  # Process each comparison
  for (i in seq_along(comparisons)) {
    comp <- comparisons[i]

    # Get result keys for both cell types
    key1 <- paste(cell_type1, comp, sep = "__")
    key2 <- paste(cell_type2, comp, sep = "__")

    # Check if both cell types have results for this comparison
    if (key1 %in% names(result_list) && key2 %in% names(result_list)) {
      # Extract features (genes or pathways) based on level
      if (level == "gene") {
        # Extract significant DEGs
        features1 <- result_list[[key1]]$gene[result_list[[key1]]$significant]
        features2 <- result_list[[key2]]$gene[result_list[[key2]]$significant]
      } else if (level == "pathway") {
        # Extract significant pathways
        tryCatch({
          features1 <- result_list[[key1]]$enrichment$ID
          features2 <- result_list[[key2]]$enrichment$ID
        }, error = function(e) {
          features1 <- character(0)
          features2 <- character(0)
        })
      }

      # Calculate similarity based on method
      if (length(features1) > 0 && length(features2) > 0) {
        comparison_similarities[i] <- calculateFeatureSimilarity(
          features1, features2, method
        )
      } else {
        comparison_similarities[i] <- 0
      }
    } else {
      # One or both cell types missing for this comparison
      comparison_similarities[i] <- 0
    }
  }

  # Return average similarity across comparisons
  if (all(is.na(comparison_similarities))) {
    return(0)
  } else {
    return(mean(comparison_similarities, na.rm = TRUE))
  }
}

#' Calculate Similarity Between Two Feature Sets
#'
#' @param features1 First feature set
#' @param features2 Second feature set
#' @param method Similarity metric
#'
#' @return Numeric similarity score
#' @keywords internal
calculateFeatureSimilarity <- function(features1,
                                       features2,
                                       method = "jaccard") {

  # Calculate similarity based on method
  if (method == "jaccard") {
    # Jaccard similarity: |A ∩ B| / |A ∪ B|
    intersection <- length(intersect(features1, features2))
    union <- length(union(features1, features2))

    if (union == 0) {
      return(0)
    } else {
      return(intersection / union)
    }

  } else if (method == "overlap") {
    # Overlap coefficient: |A ∩ B| / min(|A|, |B|)
    intersection <- length(intersect(features1, features2))
    min_size <- min(length(features1), length(features2))

    if (min_size == 0) {
      return(0)
    } else {
      return(intersection / min_size)
    }

  } else if (method == "correlation") {
    # Get all genes from both cell types
    all_features <- union(features1, features2)

    # Create binary vectors
    vec1 <- as.numeric(all_features %in% features1)
    vec2 <- as.numeric(all_features %in% features2)

    # Calculate correlation
    corr <- cor(vec1, vec2)

    if (is.na(corr)) {
      return(0)
    } else {
      return(max(0, corr))  # Ensure non-negative
    }

  } else {
    stop(paste0("Unsupported similarity method: ", method))
  }
}
