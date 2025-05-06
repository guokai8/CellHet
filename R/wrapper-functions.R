#' Wrapper function to find differentially expressed genes
#'
#' @param object Seurat or SingleCellExperiment object
#' @param group_var Name of the column in metadata containing group information
#' @param cell_type_var Name of the column in metadata containing cell type annotations
#' @param min_cells Minimum number of cells required per group for comparison
#' @param logfc.threshold Log fold-change threshold for DEGs
#' @param p_val_adj_threshold Adjusted p-value threshold for significance
#' @param p_val_cutoff p-value threshold for enrichment analysis
#' @param test.use Statistical test to use for differential expression
#' @param workers Number of cores to use for parallel processing (NULL for sequential)
#' @param future.globals.maxSize Maximum size of global variables for parallel processing
#' @param gene_id_mapping Optional data frame for gene ID mapping with columns "gene" and "id"
#' @param enrich Whether to perform enrichment analysis
#' @param annot_data Annotation database to use if enrich is TRUE
#' @param return_summary Whether to return a summary plot with the results
#' @param reference_group Optional, specify a reference group to compare all other groups against
#' @param reference_cell_type Optional, specify a reference cell type to focus analysis on
#'
#' @return A list containing DEG results for each cell type and comparison
#' @export
findDifferentialGenes <- function(object,
                                  group_var = "condition",
                                  cell_type_var = "cell_type",
                                  min_cells = 3,
                                  logfc.threshold = 0.25,
                                  p_val_adj_threshold = 0.05,
                                  p_val_cutoff = 0.05,
                                  test.use = "wilcox",
                                  workers = NULL,
                                  future.globals.maxSize = 1024^3,
                                  gene_id_mapping = NULL,
                                  enrich = FALSE,
                                  annot_data = NULL,
                                  use_gsea = FALSE,
                                  return_summary = TRUE,
                                  reference_group = NULL,
                                  reference_cell_type = NULL) {

  # Check if input is valid
  is_seurat <- inherits(object, "Seurat")
  is_sce <- inherits(object, "SingleCellExperiment")

  if (!is_seurat && !is_sce) {
    stop("Input must be a Seurat or SingleCellExperiment object")
  }

  # Extract metadata
  if (is_seurat) {
    metadata <- object@meta.data
  } else {
    metadata <- as.data.frame(colData(object))
  }

  # Check if required columns exist in metadata
  if (!group_var %in% colnames(metadata)) {
    stop(paste0("Column '", group_var, "' not found in metadata"))
  }

  if (!cell_type_var %in% colnames(metadata)) {
    stop(paste0("Column '", cell_type_var, "' not found in metadata"))
  }

  # Check reference group if specified
  if (!is.null(reference_group)) {
    if (!reference_group %in% metadata[[group_var]]) {
      stop(paste0("Reference group '", reference_group, "' not found in '", group_var, "' column"))
    }
  }

  # Check reference cell type if specified
  if (!is.null(reference_cell_type)) {
    if (!reference_cell_type %in% metadata[[cell_type_var]]) {
      stop(paste0("Reference cell type '", reference_cell_type, "' not found in '", cell_type_var, "' column"))
    }
  }

  # Set up parallelization if requested
  if (!is.null(workers)) {
    if (!requireNamespace("future", quietly = TRUE)) {
      warning("Package 'future' needed for parallel processing. Falling back to sequential processing.")
    } else {
      message(paste0("Setting up parallel processing with ", workers, " workers"))
      future::plan(future::multicore, workers = workers)
      options(future.globals.maxSize = future.globals.maxSize)
    }
  }

  # Run differential expression analysis
  message("Running differential expression analysis...")

  # Create custom comparisons if reference group is provided
  custom_comparisons <- NULL
  if (!is.null(reference_group)) {
    # Get all other groups
    other_groups <- unique(metadata[[group_var]])
    other_groups <- other_groups[other_groups != reference_group]

    # Create comparisons between reference and each other group
    custom_comparisons <- lapply(other_groups, function(g) c(g, reference_group))
    message(paste0("Using '", reference_group, "' as reference group for all comparisons"))
  }

  # Filter for specific cell type if reference_cell_type is provided
  cell_type_subset <- NULL
  if (!is.null(reference_cell_type)) {
    message(paste0("Using '", reference_cell_type, "' as reference cell type"))
    # We'll handle this after running compareDEGs
  }

  # Run differential expression analysis
  deg_results <- compareDEGs(
    sce = object,
    group_var = group_var,
    cell_type_var = cell_type_var,
    min_cells = min_cells,
    test.use = test.use,
    logfc.threshold = logfc.threshold,
    p_val_adj_threshold = p_val_adj_threshold,
    custom_comparisons = custom_comparisons,
    cores = if (!is.null(workers)) workers else 1
  )

  # Extract all DEGs into a combined data frame
  message("Combining results...")
  all_degs <- data.frame()

  for (result_key in names(deg_results$degs)) {
    # Extract cell type and comparison from key
    key_parts <- strsplit(result_key, "__")[[1]]
    cell_type <- key_parts[1]
    comparison <- key_parts[2]

    # Skip if we're only interested in reference cell type
    if (!is.null(reference_cell_type) && cell_type != reference_cell_type) {
      next
    }

    # Get DEG results
    degs <- deg_results$degs[[result_key]]

    # Add cell type and comparison columns
    degs$cell_type <- cell_type
    degs$comparison <- comparison

    # Split comparison into groups
    comp_parts <- strsplit(comparison, "_vs_")[[1]]
    degs$group1 <- comp_parts[1]
    degs$group2 <- comp_parts[2]

    # Add reference status
    if (!is.null(reference_group)) {
      degs$is_reference <- degs$group2 == reference_group
    }

    # Append to combined data frame
    all_degs <- rbind(all_degs, degs)
  }

  # Add gene ID mapping if provided
  if (!is.null(gene_id_mapping)) {
    message("Adding gene ID annotations...")
    if (all(c("gene", "id") %in% colnames(gene_id_mapping))) {
      all_degs <- merge(all_degs, gene_id_mapping, by = "gene", all.x = TRUE)
    } else {
      warning("gene_id_mapping must contain 'gene' and 'id' columns. Skipping ID mapping.")
    }
  }

  # Run pathway analysis if requested
  if (enrich) {
    message("Running pathway analysis...")
    if (!exists("gene_sets")) {
      warning("No gene sets available for pathway analysis. Skipping pathway annotation.")
    } else {

      pathway_results <- runPathwayAnalysis(
        deg_results = deg_results,
        annot_data = annot_data,
        pvalue = p_val_cutoff,
        use_gsea = use_gsea
      )

      # Add pathway information to results
      deg_results$pathway_results <- pathway_results
    }
  }

  # Update the summary to account for both significance and logFC threshold
  # Check if we have the passes_lfc_threshold column
  has_lfc_threshold_info <- "passes_lfc_threshold" %in% colnames(all_degs)

  if (has_lfc_threshold_info) {
    # Create an updated summary
    summary_data <- all_degs %>%
      dplyr::group_by(cell_type, comparison) %>%
      dplyr::summarize(
        n_up = sum(significant & direction == "up" & passes_lfc_threshold),
        n_down = sum(significant & direction == "down" & passes_lfc_threshold),
        n_cells_group1 = dplyr::first(n_cells_group1),
        n_cells_group2 = dplyr::first(n_cells_group2),
        n_genes_tested = n(),
        n_significant = sum(significant & passes_lfc_threshold),
        group1 = dplyr::first(group1),
        group2 = dplyr::first(group2),
        .groups = "drop"
      )

    deg_results$summary <- summary_data
  }

  # Create summary visualization if requested
  if (return_summary) {
    message("Generating summary visualization...")
    summary_plot <- summarizeDEGs(
      deg_results = deg_results,
      plot_type = "heatmap",
      direction = "both"
    )
    deg_results$summary_plot <- summary_plot
  }

  # Add the combined data frame to results
  deg_results$all_degs <- all_degs

  # Return results
  message("Done!")
  return(deg_results)
}

#' Quick Comparison of Two Specific Groups for a Cell Type
#'
#' A simplified wrapper to quickly compare two specific groups for a given cell type
#' or set of cell types.
#'
#' @param object Seurat or SingleCellExperiment object
#' @param cell_types Character vector of cell types to analyze (or "all" for all cell types)
#' @param group1 Name of first group to compare
#' @param group2 Name of second group to compare
#' @param group_var Name of the column in metadata containing group information
#' @param cell_type_var Name of the column in metadata containing cell type annotations
#' @param logfc.threshold Log fold-change threshold for DEG identification
#' @param p_val_adj_threshold Adjusted p-value threshold for significance
#' @param test.use Statistical test to use for differential expression
#' @param min_cells Minimum number of cells required per group for comparison
#' @param gene_id_mapping Optional data frame for gene ID mapping
#' @param return_plot Whether to return a volcano plot of the results
#' @param reference_group Optional, specify which group should be used as reference
#'
#' @return A data frame of differential genes or a list with results and plot
#' @export
quickCompare <- function(object,
                         cell_types = "all",
                         group1,
                         group2,
                         group_var = "condition",
                         cell_type_var = "cell_type",
                         logfc.threshold = 0.25,
                         p_val_adj_threshold = 0.05,
                         test.use = "wilcox",
                         min_cells = 3,
                         gene_id_mapping = NULL,
                         return_plot = TRUE,
                         reference_group = NULL) {

  # Check if input is valid
  is_seurat <- inherits(object, "Seurat")
  is_sce <- inherits(object, "SingleCellExperiment")

  if (!is_seurat && !is_sce) {
    stop("Input must be a Seurat or SingleCellExperiment object")
  }

  # Extract metadata
  if (is_seurat) {
    metadata <- object@meta.data
  } else {
    metadata <- as.data.frame(colData(object))
  }

  # Check if required columns exist in metadata
  if (!group_var %in% colnames(metadata)) {
    stop(paste0("Column '", group_var, "' not found in metadata"))
  }

  if (!cell_type_var %in% colnames(metadata)) {
    stop(paste0("Column '", cell_type_var, "' not found in metadata"))
  }

  # Check if groups exist
  if (!group1 %in% metadata[[group_var]]) {
    stop(paste0("Group '", group1, "' not found in '", group_var, "' column"))
  }

  if (!group2 %in% metadata[[group_var]]) {
    stop(paste0("Group '", group2, "' not found in '", group_var, "' column"))
  }

  # Get cell types to analyze
  available_cell_types <- unique(metadata[[cell_type_var]])

  if (identical(cell_types, "all")) {
    cell_types <- available_cell_types
  } else {
    # Check if specified cell types exist
    invalid_cell_types <- cell_types[!cell_types %in% available_cell_types]
    if (length(invalid_cell_types) > 0) {
      warning(paste0("The following cell types were not found: ",
                     paste(invalid_cell_types, collapse = ", "),
                     ". They will be skipped."))
      cell_types <- cell_types[cell_types %in% available_cell_types]
    }

    if (length(cell_types) == 0) {
      stop("No valid cell types specified")
    }
  }

  # Define custom comparison
  custom_comparison <- list(c(group1, group2))

  # Results storage
  all_results <- list()

  # Process each cell type
  message(paste0("Comparing ", group1, " vs ", group2, " for ", length(cell_types), " cell types..."))

  for (ct in cell_types) {
    message(paste0("Processing ", ct, "..."))

    # Subset data for current cell type
    if (is_seurat) {
      cell_subset <- subset(object, cells = which(metadata[[cell_type_var]] == ct))

      # Check cell counts
      n_cells_group1 <- sum(cell_subset@meta.data[[group_var]] == group1)
      n_cells_group2 <- sum(cell_subset@meta.data[[group_var]] == group2)
    } else {
      cell_subset <- object[, metadata[[cell_type_var]] == ct]

      # Check cell counts
      n_cells_group1 <- sum(colData(cell_subset)[[group_var]] == group1)
      n_cells_group2 <- sum(colData(cell_subset)[[group_var]] == group2)
    }

    # Skip if not enough cells
    if (n_cells_group1 < min_cells || n_cells_group2 < min_cells) {
      message(paste0("Skipping ", ct, ": insufficient cells (",
                     group1, ": ", n_cells_group1, ", ",
                     group2, ": ", n_cells_group2, ")"))
      next
    }

    # Run differential expression
    if (is_seurat) {
      # Set identity to group variable
      Seurat::Idents(cell_subset) <- cell_subset@meta.data[[group_var]]

      # Determine which group is reference (ident.2)
      ident1 <- group1
      ident2 <- group2

      # Override if reference_group is specified
      if (!is.null(reference_group)) {
        if (reference_group == group1) {
          ident1 <- group2
          ident2 <- group1
        } else if (reference_group == group2) {
          # Default behavior already correct
        } else {
          warning(paste0("Reference group '", reference_group,
                         "' is not part of the comparison. Using default: ", group2, " as reference."))
        }
      }

      # [REVISED] Run FindMarkers with no logFC threshold to get all genes
      de_results <- Seurat::FindMarkers(
        object = cell_subset,
        ident.1 = ident1,
        ident.2 = ident2,
        test.use = test.use,
        logfc.threshold = 0,  # [REVISED] No threshold to get all genes
        min.pct = 0  # [REVISED] Include all genes
      )

      # [REVISED] Mark significant genes and those passing logFC threshold
      de_results$significant <- de_results$p_val_adj < p_val_adj_threshold
      de_results$passes_lfc_threshold <- abs(de_results$avg_log2FC) >= logfc.threshold

    } else {
      # For SingleCellExperiment, use scran
      if (!requireNamespace("scran", quietly = TRUE)) {
        stop("Package 'scran' needed for SingleCellExperiment objects. Please install it.")
      }

      # Create factor for comparison
      group_factor <- factor(colData(cell_subset)[[group_var]])
      group_factor <- droplevels(group_factor)

      # [REVISED] Run findMarkers with no logFC threshold
      design <- model.matrix(~group_factor)
      de_results <- scran::findMarkers(
        cell_subset,
        groups = group_factor,
        design = design,
        direction = "any",
        lfc = 0  # [REVISED] No threshold to get all genes
      )

      # Format results to match Seurat output
      de_results <- as.data.frame(de_results[[group1]])
      de_results$p_val_adj <- de_results$FDR
      de_results$avg_log2FC <- de_results$logFC

      # [REVISED] Mark significant genes and those passing logFC threshold
      de_results$significant <- de_results$p_val_adj < p_val_adj_threshold
      de_results$passes_lfc_threshold <- abs(de_results$avg_log2FC) >= logfc.threshold
    }

    # Add gene column
    de_results$gene <- rownames(de_results)

    # Add direction
    de_results$direction <- ifelse(de_results$avg_log2FC > 0, "up", "down")

    # Add cell type information
    de_results$cell_type <- ct
    de_results$comparison <- paste(group1, "vs", group2, sep = "_")

    # Add gene ID mapping if provided
    if (!is.null(gene_id_mapping)) {
      if (all(c("gene", "id") %in% colnames(gene_id_mapping))) {
        de_results <- merge(de_results, gene_id_mapping, by = "gene", all.x = TRUE)
      }
    }

    # Store results
    all_results[[ct]] <- de_results
  }

  # Combine all results
  if (length(all_results) == 0) {
    stop("No valid comparisons could be performed")
  }

  combined_results <- do.call(rbind, all_results)

  # Create volcano plot if requested
  if (return_plot && requireNamespace("ggplot2", quietly = TRUE)) {
    plot_list <- list()

    for (ct in names(all_results)) {
      # Get data for this cell type
      plot_data <- all_results[[ct]]

      # [REVISED] Create labels for top genes - use genes that pass both significance and logFC threshold
      plot_data$label <- ""
      if (nrow(plot_data) > 0) {
        # For volcano plot, highlight genes that pass both thresholds
        is_deg <- plot_data$significant & plot_data$passes_lfc_threshold

        # Label top up and down genes
        top_up <- plot_data$gene[is_deg & plot_data$avg_log2FC > 0]
        top_up <- top_up[order(plot_data$p_val_adj[is_deg & plot_data$avg_log2FC > 0])]
        top_up <- head(top_up, 10)

        top_down <- plot_data$gene[is_deg & plot_data$avg_log2FC < 0]
        top_down <- top_down[order(plot_data$p_val_adj[is_deg & plot_data$avg_log2FC < 0])]
        top_down <- head(top_down, 10)

        plot_data$label[plot_data$gene %in% c(top_up, top_down)] <- plot_data$gene[plot_data$gene %in% c(top_up, top_down)]
      }

      # [REVISED] Create volcano plot - highlight genes that pass both significance and logFC threshold
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = avg_log2FC, y = -log10(p_val_adj),
                                                   color = direction, label = label)) +
        ggplot2::geom_point(ggplot2::aes(alpha = significant & passes_lfc_threshold)) + # [REVISED]
        ggplot2::geom_vline(xintercept = c(-logfc.threshold, logfc.threshold),
                            linetype = "dashed", color = "gray50") +
        ggplot2::geom_hline(yintercept = -log10(p_val_adj_threshold),
                            linetype = "dashed", color = "gray50") +
        ggplot2::scale_color_manual(values = c("down" = "deepskyblue", "up" = "darkorange")) +
        ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
        ggplot2::labs(
          title = paste0(ct, ": ", group1, " vs ", group2),
          subtitle = paste0(sum(plot_data$significant & plot_data$passes_lfc_threshold), # [REVISED]
                            " significant DEGs (",
                            sum(plot_data$significant & plot_data$passes_lfc_threshold & # [REVISED]
                                  plot_data$direction == "up"), " up, ",
                            sum(plot_data$significant & plot_data$passes_lfc_threshold & # [REVISED]
                                  plot_data$direction == "down"), " down)"),
          x = "Log2 Fold Change",
          y = "-log10(adjusted p-value)"
        ) +
        ggplot2::theme_bw()

      # Add labels if there are significant genes
      if (sum(plot_data$significant & plot_data$passes_lfc_threshold) > 0) { # [REVISED]
        if (requireNamespace("ggrepel", quietly = TRUE)) {
          p <- p + ggrepel::geom_text_repel(
            data = subset(plot_data, label != ""),
            box.padding = 0.5,
            max.overlaps = 20
          )
        } else {
          # Fallback if ggrepel not available
          p <- p + ggplot2::geom_text(data = subset(plot_data, label != ""), hjust = 0, vjust = 0)
        }
      }

      plot_list[[ct]] <- p
    }

    # Return results and plots
    return(list(
      results = combined_results,
      plots = plot_list
    ))
  } else {
    # Return just the results
    return(combined_results)
  }
}

#' Export Differential Gene Results to Excel
#'
#' @param deg_results Output from findDifferentialGenes or quickCompare
#' @param file_path Path to save the Excel file
#' @param split_by How to split results into worksheets ("cell_type", "comparison", or "none")
#' @param include_pathways Whether to include pathway analysis results if available
#'
#' @return Invisibly returns the file path
#' @export
exportDEGResults <- function(deg_results,
                             file_path,
                             split_by = "cell_type",
                             include_pathways = TRUE) {

  # Check for required package
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' needed for Excel export. Please install it.")
  }

  # Extract results
  if (is.data.frame(deg_results)) {
    # Direct data frame input (from quickCompare)
    all_degs <- deg_results
  } else if ("all_degs" %in% names(deg_results)) {
    # Output from findDifferentialGenes
    all_degs <- deg_results$all_degs
  } else if ("degs" %in% names(deg_results)) {
    # Output from compareDEGs
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
    stop("Unrecognized input format. Must be a data frame or output from findDifferentialGenes/compareDEGs.")
  }

  # Create a new workbook
  wb <- openxlsx::createWorkbook()

  # Create a summary sheet
  openxlsx::addWorksheet(wb, "Summary")

  # Check for the passes_lfc_threshold column
  has_lfc_threshold_info <- "passes_lfc_threshold" %in% colnames(all_degs)

  # Add summary table
  if (has_lfc_threshold_info) {
    summary_table <- all_degs %>%
      dplyr::group_by(cell_type, comparison) %>%
      dplyr::summarize(
        total_genes = n(),
        significant_genes = sum(significant & passes_lfc_threshold), # [REVISED]
        up_regulated = sum(significant & direction == "up" & passes_lfc_threshold), # [REVISED]
        down_regulated = sum(significant & direction == "down" & passes_lfc_threshold), # [REVISED]
        min_logFC = min(avg_log2FC, na.rm = TRUE),
        max_logFC = max(avg_log2FC, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    summary_table <- all_degs %>%
      dplyr::group_by(cell_type, comparison) %>%
      dplyr::summarize(
        total_genes = n(),
        significant_genes = sum(significant),
        up_regulated = sum(significant & direction == "up"),
        down_regulated = sum(significant & direction == "down"),
        min_logFC = min(avg_log2FC, na.rm = TRUE),
        max_logFC = max(avg_log2FC, na.rm = TRUE),
        .groups = "drop"
      )
  }

  openxlsx::writeData(wb, "Summary", summary_table)

  # Create styles
  headerStyle <- openxlsx::createStyle(
    fontColour = "#FFFFFF", fgFill = "#4F81BD",
    halign = "center", valign = "center", textDecoration = "bold"
  )

  # Apply styles to summary sheet
  openxlsx::addStyle(wb, "Summary", headerStyle, rows = 1, cols = 1:ncol(summary_table))

  # Create sheets based on split_by parameter
  if (split_by == "cell_type") {
    # Split by cell type
    for (ct in unique(all_degs$cell_type)) {
      # Create sheet name (sanitize if needed)
      sheet_name <- gsub("[\\/:*?\"<>|]", "_", ct)

      # Add worksheet
      openxlsx::addWorksheet(wb, sheet_name)

      # Filter data for this cell type
      ct_data <- all_degs[all_degs$cell_type == ct, ]

      # Write data
      openxlsx::writeData(wb, sheet_name, ct_data)

      # Apply styles
      openxlsx::addStyle(wb, sheet_name, headerStyle, rows = 1, cols = 1:ncol(ct_data))
    }
  } else if (split_by == "comparison") {
    # Split by comparison
    for (comp in unique(all_degs$comparison)) {
      # Create sheet name (sanitize if needed)
      sheet_name <- gsub("[\\/:*?\"<>|]", "_", comp)

      # Add worksheet
      openxlsx::addWorksheet(wb, sheet_name)

      # Filter data for this comparison
      comp_data <- all_degs[all_degs$comparison == comp, ]

      # Write data
      openxlsx::writeData(wb, sheet_name, comp_data)

      # Apply styles
      openxlsx::addStyle(wb, sheet_name, headerStyle, rows = 1, cols = 1:ncol(comp_data))
    }
  } else {
    # No splitting, just one sheet with all results
    openxlsx::addWorksheet(wb, "All_DEGs")
    openxlsx::writeData(wb, "All_DEGs", all_degs)
    openxlsx::addStyle(wb, "All_DEGs", headerStyle, rows = 1, cols = 1:ncol(all_degs))
  }

  # Include pathway results if available and requested
  if (include_pathways &&
      "pathway_results" %in% names(deg_results) &&
      "pathway_results" %in% names(deg_results$pathway_results)) {

    # Add pathway summary sheet
    openxlsx::addWorksheet(wb, "Pathway_Summary")

    # Extract pathway results
    pathway_list <- deg_results$pathway_results$pathway_results

    # Create summary table
    pathway_summary <- data.frame()

    for (result_key in names(pathway_list)) {
      # Extract cell type and comparison from key
      key_parts <- strsplit(result_key, "__")[[1]]
      cell_type <- key_parts[1]
      comparison <- key_parts[2]

      # Get pathway results
      pathway_result <- pathway_list[[result_key]]

      # Extract enrichment results
      enrichment <- pathway_result$enrichment

      # Count pathways
      n_pathways <- nrow(enrichment)

      # Add to summary
      if (n_pathways > 0) {
        pathway_summary <- rbind(pathway_summary, data.frame(
          cell_type = cell_type,
          comparison = comparison,
          n_pathways = n_pathways,
          top_pathway = ifelse(n_pathways > 0, enrichment$Description[1], NA),
          top_pvalue = ifelse(n_pathways > 0, enrichment$pvalue[1], NA)
        ))
      }
    }

    # Write pathway summary
    if (nrow(pathway_summary) > 0) {
      openxlsx::writeData(wb, "Pathway_Summary", pathway_summary)
      openxlsx::addStyle(wb, "Pathway_Summary", headerStyle, rows = 1, cols = 1:ncol(pathway_summary))

      # Add detailed pathway sheets
      for (result_key in names(pathway_list)) {
        # Extract cell type and comparison from key
        key_parts <- strsplit(result_key, "__")[[1]]
        cell_type <- key_parts[1]
        comparison <- key_parts[2]

        # Create sheet name
        sheet_name <- paste0("Pathway_", gsub("[\\/:*?\"<>|]", "_", cell_type), "_",
                             gsub("[\\/:*?\"<>|]", "_", comparison))

        # Limit sheet name length
        if (nchar(sheet_name) > 31) {
          sheet_name <- substr(sheet_name, 1, 31)
        }

        # Get pathway results
        pathway_result <- pathway_list[[result_key]]

        # Extract enrichment results
        enrichment <- pathway_result$enrichment

        # Add worksheet if pathways exist
        if (nrow(enrichment) > 0) {
          # Add worksheet
          openxlsx::addWorksheet(wb, sheet_name)

          # Write data
          openxlsx::writeData(wb, sheet_name, enrichment)

          # Apply styles
          openxlsx::addStyle(wb, sheet_name, headerStyle, rows = 1, cols = 1:ncol(enrichment))
        }
      }
    }
  }

  # Save workbook
  openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)

  message(paste0("Results exported to ", file_path))

  # Return file path invisibly
  invisible(file_path)
}
