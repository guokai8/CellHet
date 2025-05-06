#' Create a Dotplot Visualization of Pathway Enrichment Results from runPathwayAnalysis
#'
#' @param pathway_results Output from runPathwayAnalysis function
#' @param top_n Number of top pathways to include (default: 10)
#' @param color_by Variable to color dots by - "RichFactor", "Padj", "Pvalue", etc. (default: "RichFactor")
#' @param size_by Variable to size dots by - "Pvalue", "Padj", or "Significant" (default: "Pvalue")
#' @param show_comparison_names Whether to show full comparison names or simplify them (default: TRUE)
#' @param flip Flip the coordinates to have pathways on the y-axis (default: TRUE)
#' @param title Custom title for the plot (default: "Pathway Enrichment Analysis")
#' @param show_category How many categories to show (default: 10)
#'
#' @return A ggplot object visualizing pathway enrichment as a dotplot
#' @export
visualizePathwaysDotplot <- function(pathway_results,
                                     top_n = 10,
                                     color_by = "RichFactor",
                                     size_by = "Pvalue",
                                     show_comparison_names = TRUE,
                                     flip = TRUE,
                                     title = "Pathway Enrichment Analysis",
                                     show_category = 10) {

  # Check if ggplot2 is installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for visualizing pathways. Please install it.")
  }

  # Extract pathway results from runPathwayAnalysis output
  if (!inherits(pathway_results, "list") || !"pathway_results" %in% names(pathway_results)) {
    stop("Input must be output from runPathwayAnalysis function")
  }

  # Initialize data frame for combined results
  all_results <- data.frame()

  # Extract pathway results
  pathway_list <- pathway_results$pathway_results

  if (length(pathway_list) == 0) {
    stop("No pathway results found in the input data")
  }

  # Extract enrichment results into a combined data frame
  for (comp_name in names(pathway_list)) {
    # Get pathway results
    pathway_result <- pathway_list[[comp_name]]

    # Extract enrichment results
    enrichment <- pathway_result$enrichment

    # Skip if no results
    if (nrow(enrichment) == 0) {
      next
    }

    # Add comparison information
    enrichment$comparison <- comp_name

    # Add to combined results
    all_results <- rbind(all_results, enrichment)
  }

  # Check if we have any results
  if (nrow(all_results) == 0) {
    stop("No significant pathway enrichment found")
  }

  # For each comparison, get top_n pathways by p-value
  top_results <- data.frame()

  for (comp in unique(all_results$comparison)) {
    comp_results <- all_results[all_results$comparison == comp, ]
    comp_results <- comp_results[order(comp_results$Pvalue), ]

    # Take top_n or all if fewer
    n_pathways <- min(show_category, nrow(comp_results))
    top_results <- rbind(top_results, comp_results[1:n_pathways, ])
  }

  # Extract values for dot size and color

  # Prepare size variable
  if (size_by == "Pvalue" && "Pvalue" %in% colnames(top_results)) {
    top_results$size_var <- -log10(top_results$Pvalue)
    size_label <- "-log10(p-value)"
  } else if (size_by == "Padj" && "Padj" %in% colnames(top_results)) {
    top_results$size_var <- -log10(top_results$Padj)
    size_label <- "-log10(adjusted p-value)"
  } else if (size_by == "Significant" && "Significant" %in% colnames(top_results)) {
    top_results$size_var <- top_results$Significant
    size_label <- "Gene Count"
  } else {
    # Default to Pvalue if the requested column doesn't exist
    if ("Pvalue" %in% colnames(top_results)) {
      top_results$size_var <- -log10(top_results$Pvalue)
      size_label <- "-log10(p-value)"
    } else if ("pvalue" %in% colnames(top_results)) {
      top_results$size_var <- -log10(top_results$pvalue)
      size_label <- "-log10(p-value)"
    } else {
      stop("Could not find p-value column in enrichment results")
    }
    warning("Requested size_by column not found, using p-value instead")
  }

  # Prepare color variable
  if (color_by == "RichFactor" && "RichFactor" %in% colnames(top_results)) {
    top_results$color_var <- top_results$RichFactor
    color_label <- "Rich Factor"
  } else if (color_by == "Padj" && "Padj" %in% colnames(top_results)) {
    top_results$color_var <- -log10(top_results$Padj)
    color_label <- "-log10(adjusted p-value)"
  } else if (color_by == "FoldEnrichment" && "FoldEnrichment" %in% colnames(top_results)) {
    top_results$color_var <- top_results$FoldEnrichment
    color_label <- "Fold Enrichment"
  } else {
    # Default to RichFactor if available, otherwise use FoldEnrichment
    if ("RichFactor" %in% colnames(top_results)) {
      top_results$color_var <- top_results$RichFactor
      color_label <- "Rich Factor"
    } else if ("FoldEnrichment" %in% colnames(top_results)) {
      top_results$color_var <- top_results$FoldEnrichment
      color_label <- "Fold Enrichment"
    } else {
      # If neither exists, use Pvalue
      top_results$color_var <- -log10(top_results$Pvalue)
      color_label <- "-log10(p-value)"
    }
    warning("Requested color_by column not found, using appropriate alternative")
  }

  # Create a more readable pathway label
  if ("Term" %in% colnames(top_results)) {
    # For richKEGG format
    top_results$pathway_label <- top_results$Term
  } else if ("Description" %in% colnames(top_results)) {
    # For GO format
    top_results$pathway_label <- top_results$Description
  } else {
    # Fallback to ID
    top_results$pathway_label <- top_results$ID
  }

  # Ensure pathway labels are factors in the right order
  if ("Pvalue" %in% colnames(top_results)) {
    pathway_order <- unique(top_results$pathway_label[order(top_results$Pvalue)])
  } else {
    pathway_order <- unique(top_results$pathway_label[order(top_results$pvalue)])
  }

  if (flip) {
    # For horizontal bars (pathways on y-axis)
    top_results$pathway_label <- factor(top_results$pathway_label, levels = rev(pathway_order))
  } else {
    # For vertical bars (pathways on x-axis)
    top_results$pathway_label <- factor(top_results$pathway_label, levels = pathway_order)
  }

  # Extract cell type and comparison if they exist in the comparison column
  if ("comparison" %in% colnames(top_results) &&
      any(grepl("__", top_results$comparison))) {

    split_comparisons <- strsplit(top_results$comparison, "__")
    top_results$cell_type <- sapply(split_comparisons, function(x) x[1])
    top_results$comp_name <- sapply(split_comparisons, function(x) x[2])

    # Use facets if both cell type and comparison exist and vary
    if (length(unique(top_results$cell_type)) > 1 && length(unique(top_results$comp_name)) > 1) {
      use_facets <- TRUE
    } else {
      use_facets <- FALSE
    }

    # Simplify comparison names if requested
    if (!show_comparison_names) {
      if (use_facets) {
        # If using facets, use only the comp_name for the x-axis
        top_results$plot_comparison <- top_results$comp_name
      } else {
        # If cell types vary but comparisons don't, use cell types
        if (length(unique(top_results$cell_type)) > 1 && length(unique(top_results$comp_name)) == 1) {
          top_results$plot_comparison <- top_results$cell_type
        }
        # If comparisons vary but cell types don't, use comparisons
        else if (length(unique(top_results$cell_type)) == 1 && length(unique(top_results$comp_name)) > 1) {
          top_results$plot_comparison <- top_results$comp_name
        }
        # If both vary but we're not using facets, use the full comparison
        else {
          top_results$plot_comparison <- top_results$comparison
        }
      }
    } else {
      # Use full comparison names
      top_results$plot_comparison <- top_results$comparison
    }
  } else {
    use_facets <- FALSE
    top_results$plot_comparison <- top_results$comparison
  }

  # Create dotplot
  library(ggplot2)

  if (use_facets) {
    # Plot with facets for multiple cell types
    if (flip) {
      p <- ggplot(top_results,
                  aes(x = comp_name, y = pathway_label,
                      size = size_var,
                      color = color_var)) +
        geom_point() +
        facet_wrap(~ cell_type, scales = "free_y") +
        scale_size_continuous(name = size_label) +
        scale_color_viridis_c(name = color_label) +
        labs(title = title,
             x = "Comparison",
             y = "Pathway") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              strip.background = element_rect(fill = "lightblue", color = "black"),
              strip.text = element_text(face = "bold"))
    } else {
      p <- ggplot(top_results,
                  aes(x = pathway_label, y = comp_name,
                      size = size_var,
                      color = color_var)) +
        geom_point() +
        facet_wrap(~ cell_type, scales = "free_x") +
        scale_size_continuous(name = size_label) +
        scale_color_viridis_c(name = color_label) +
        labs(title = title,
             y = "Comparison",
             x = "Pathway") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              strip.background = element_rect(fill = "lightblue", color = "black"),
              strip.text = element_text(face = "bold"))
    }
  } else {
    # Standard plot for single comparison or multiple comparisons without faceting
    if (flip) {
      p <- ggplot(top_results,
                  aes(x = plot_comparison, y = pathway_label,
                      size = size_var,
                      color = color_var)) +
        geom_point() +
        scale_size_continuous(name = size_label) +
        scale_color_viridis_c(name = color_label) +
        labs(title = title,
             x = "Comparison",
             y = "Pathway") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    } else {
      p <- ggplot(top_results,
                  aes(x = pathway_label, y = plot_comparison,
                      size = size_var,
                      color = color_var)) +
        geom_point() +
        scale_size_continuous(name = size_label) +
        scale_color_viridis_c(name = color_label) +
        labs(title = title,
             y = "Comparison",
             x = "Pathway") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  }

  return(p)
}

#' Save Pathway Dotplot to a File
#'
#' @param plot ggplot object from visualizePathwaysDotplot
#' @param file Output file path (default: "pathway_dotplot.pdf")
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @param dpi Resolution for raster outputs (default: 300)
#'
#' @return Invisibly returns the plot object
#' @export
savePathwayDotplot <- function(plot,
                               file = "pathway_dotplot.pdf",
                               width = 10,
                               height = 8,
                               dpi = 300) {

  # Determine file type from extension
  file_ext <- tolower(tools::file_ext(file))

  # Save file based on type
  if (file_ext == "pdf") {
    pdf(file, width = width, height = height)
    print(plot)
    dev.off()
  } else if (file_ext %in% c("png", "jpg", "jpeg", "tiff")) {
    ggsave(filename = file,
           plot = plot,
           width = width,
           height = height,
           dpi = dpi)
  } else {
    warning("Unrecognized file extension. Saving as PDF.")
    pdf("pathway_dotplot.pdf", width = width, height = height)
    print(plot)
    dev.off()
  }

  # Return plot invisibly
  invisible(plot)
}
