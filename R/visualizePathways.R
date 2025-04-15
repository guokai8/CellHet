# visualizePathways.R

#' Visualize Pathway Enrichment Results
#'
#' @param pathway_results Output from runPathwayAnalysis function
#' @param plot_type Type of visualization ("dotplot", "heatmap", "network")
#' @param top_n Number of top pathways to include per comparison
#' @param common_only If TRUE, only show pathways found in multiple comparisons
#'
#' @return A plot object visualizing pathway enrichment
#' @export
visualizePathways <- function(pathway_results,
                              plot_type = "dotplot",
                              top_n = 10,
                              common_only = FALSE) {

  # Extract pathway results
  pathway_list <- pathway_results$pathway_results

  if (length(pathway_list) == 0) {
    stop("No pathway results found in the input data")
  }

  # Extract enrichment results into a combined data frame
  all_results <- data.frame()

  for (comparison_name in names(pathway_list)) {
    # Get pathway results
    pathway_result <- pathway_list[[comparison_name]]

    # Extract enrichment results
    enrichment <- pathway_result$enrichment

    # Skip if no results
    if (nrow(enrichment) == 0) {
      next
    }

    # Add comparison information
    enrichment$comparison <- comparison_name

    # Add to combined results
    all_results <- rbind(all_results, enrichment)
  }

  # Check if we have any results
  if (nrow(all_results) == 0) {
    stop("No significant pathway enrichment found")
  }

  # If common_only, filter to pathways found in multiple comparisons
  if (common_only && length(unique(all_results$comparison)) > 1) {
    # Count occurrences of each pathway
    pathway_counts <- table(all_results$ID)

    # Keep only pathways found in multiple comparisons
    common_pathways <- names(pathway_counts[pathway_counts > 1])

    # Filter results
    all_results <- all_results[all_results$ID %in% common_pathways, ]

    if (nrow(all_results) == 0) {
      stop("No common pathways found across comparisons")
    }
  }

  # For each comparison, get top_n pathways by p-value
  top_results <- data.frame()

  for (comp in unique(all_results$comparison)) {
    comp_results <- all_results[all_results$comparison == comp, ]
    comp_results <- comp_results[order(comp_results$pvalue), ]

    # Take top_n or all if fewer
    n_pathways <- min(top_n, nrow(comp_results))
    top_results <- rbind(top_results, comp_results[1:n_pathways, ])
  }

  # Create appropriate visualization
  if (plot_type == "dotplot") {
    # Create dotplot
    # Extract -log10(p-value) for dot size
    top_results$neg_log_p <- -log10(top_results$pvalue)

    # Extract gene ratio for dot color
    top_results$gene_ratio <- sapply(strsplit(top_results$GeneRatio, "/"),
                                     function(x) as.numeric(x[1]) / as.numeric(x[2]))

    # Create a more readable pathway label
    top_results$pathway_label <- top_results$Description

    # Ensure pathway labels are factors in the right order
    pathway_order <- unique(top_results$pathway_label[order(top_results$pvalue)])
    top_results$pathway_label <- factor(top_results$pathway_label, levels = rev(pathway_order))

    p <- ggplot(top_results,
                aes(x = comparison, y = pathway_label,
                    size = neg_log_p,
                    color = gene_ratio)) +
      geom_point() +
      scale_size_continuous(name = "-log10(p-value)") +
      scale_color_viridis_c(name = "Gene Ratio") +
      labs(title = "Pathway Enrichment Analysis",
           x = "Comparison",
           y = "Pathway") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

  } else if (plot_type == "heatmap") {
    # Create heatmap of pathway enrichment
    # Create a wide format matrix for heatmap
    # Use -log10(p-value) as the value
    top_results$neg_log_p <- -log10(top_results$pvalue)

    # Reshape data for heatmap
    wide_data <- reshape2::dcast(top_results,
                                 Description ~ comparison,
                                 value.var = "neg_log_p",
                                 fill = 0)

    # Convert to matrix
    rownames(wide_data) <- wide_data$Description
    wide_data$Description <- NULL
    heatmap_matrix <- as.matrix(wide_data)

    if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      # Use ComplexHeatmap
      p <- ComplexHeatmap::Heatmap(
        matrix = heatmap_matrix,
        name = "-log10(p-value)",
        row_title = "Pathway",
        column_title = "Comparison",
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 10),
        column_names_rot = 45,
        col = viridis::viridis(100)
      )
    } else {
      # Fallback to ggplot2
      heatmap_data <- reshape2::melt(heatmap_matrix,
                                     varnames = c("Pathway", "Comparison"),
                                     value.name = "neg_log_p")

      p <- ggplot(heatmap_data, aes(x = Comparison, y = Pathway, fill = neg_log_p)) +
        geom_tile() +
        scale_fill_viridis_c(name = "-log10(p-value)") +
        labs(title = "Pathway Enrichment Heatmap") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }

  } else if (plot_type == "network") {
    if (!requireNamespace("igraph", quietly = TRUE) ||
        !requireNamespace("ggraph", quietly = TRUE)) {
      stop("Packages 'igraph' and 'ggraph' needed for network plots. Please install them.")
    }

    # Create a pathway-comparison network
    # Extract edges (pathway-comparison pairs)
    edges <- top_results[, c("Description", "comparison", "neg_log_p", "gene_ratio")]
    names(edges) <- c("pathway", "comparison", "neg_log_p", "gene_ratio")

    # Create nodes for pathways and comparisons
    pathway_nodes <- data.frame(
      name = unique(edges$pathway),
      type = "pathway",
      stringsAsFactors = FALSE
    )

    comparison_nodes <- data.frame(
      name = unique(edges$comparison),
      type = "comparison",
      stringsAsFactors = FALSE
    )

    # Combine nodes
    nodes <- rbind(pathway_nodes, comparison_nodes)

    # Create a graph
    g <- igraph::graph_from_data_frame(
      d = edges[, c("pathway", "comparison", "neg_log_p", "gene_ratio")],
      vertices = nodes,
      directed = FALSE
    )

    # Create a layout
    layout <- igraph::layout_with_fr(g)

    # Create the network plot
    p <- ggraph::ggraph(g, layout = "fr") +
      ggraph::geom_edge_link(ggplot2::aes(width = neg_log_p, alpha = gene_ratio),
                             edge_colour = "skyblue") +
      ggraph::geom_node_point(ggplot2::aes(color = type, size = type)) +
      ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE) +
      ggraph::scale_edge_width(range = c(0.5, 3)) +
      ggraph::scale_edge_alpha(range = c(0.3, 1)) +
      ggplot2::scale_color_manual(values = c("pathway" = "orange", "comparison" = "steelblue")) +
      ggplot2::scale_size_manual(values = c("pathway" = 5, "comparison" = 7)) +
      ggraph::theme_graph() +
      ggplot2::labs(title = "Pathway-Comparison Network",
                    edge_width = "-log10(p-value)",
                    edge_alpha = "Gene Ratio")

  } else {
    stop(paste0("Unsupported plot type: ", plot_type))
  }

  return(p)
}
