#' Compare DEGs Across Cell Types and Conditions
#'
#' @param sce SingleCellExperiment or Seurat object with annotated cell types
#' @param group_var Name of the column in metadata containing group information
#' @param cell_type_var Name of the column in metadata containing cell type annotations
#' @param min_cells Minimum number of cells required per group for comparison
#' @param test.use Statistical test to use for differential expression
#' @param logfc.threshold Log fold-change threshold for DEGs identification (not for filtering results)
#' @param p_val_adj_threshold Adjusted p-value threshold for significance
#' @param custom_comparisons Optional list of custom comparisons to perform
#' @param cores Number of cores to use for parallel processing
#'
#' @return A list containing DEG results for each cell type and comparison
#' @export
compareDEGs <- function(sce,
                        group_var = "condition",
                        cell_type_var = "cell_type",
                        min_cells = 3,
                        test.use = "wilcox",
                        logfc.threshold = 0.25,
                        p_val_adj_threshold = 0.05,
                        custom_comparisons = NULL,
                        cores = 1) {

  # Check if input is Seurat or SingleCellExperiment
  is_seurat <- inherits(sce, "Seurat")
  is_sce <- inherits(sce, "SingleCellExperiment")

  if (!is_seurat && !is_sce) {
    stop("Input must be a Seurat or SingleCellExperiment object")
  }

  # Extract metadata
  if (is_seurat) {
    metadata <- sce@meta.data
  } else {
    metadata <- as.data.frame(colData(sce))
  }

  # Check if required columns exist in metadata
  if (!group_var %in% colnames(metadata)) {
    stop(paste0("Column '", group_var, "' not found in metadata"))
  }

  if (!cell_type_var %in% colnames(metadata)) {
    stop(paste0("Column '", cell_type_var, "' not found in metadata"))
  }

  # Extract unique cell types and groups
  unique_cell_types <- unique(metadata[[cell_type_var]])
  unique_groups <- unique(metadata[[group_var]])

  # Define comparisons
  if (is.null(custom_comparisons)) {
    # Generate all pairwise comparisons
    all_comparisons <- combn(unique_groups, 2, simplify = FALSE)
  } else {
    # Use user-defined comparisons
    all_comparisons <- custom_comparisons
    # Validate custom comparisons
    for (comp in all_comparisons) {
      if (!all(comp %in% unique_groups)) {
        stop("Custom comparison contains groups not present in the data")
      }
    }
  }

  # Initialize results storage
  deg_results <- list()
  summary_stats <- data.frame()

  # Process each cell type
  for (cell_type in unique_cell_types) {
    # Subset data for current cell type
    if (is_seurat) {
      cell_subset <- subset(sce, cells = which(metadata[[cell_type_var]] == cell_type))
    } else {
      cell_subset <- sce[, metadata[[cell_type_var]] == cell_type]
    }

    # Process each comparison
    for (comp_idx in seq_along(all_comparisons)) {
      comparison <- all_comparisons[[comp_idx]]
      group1 <- comparison[1]
      group2 <- comparison[2]

      # Check if enough cells for comparison
      if (is_seurat) {
        n_cells_group1 <- sum(cell_subset@meta.data[[group_var]] == group1)
        n_cells_group2 <- sum(cell_subset@meta.data[[group_var]] == group2)
      } else {
        n_cells_group1 <- sum(colData(cell_subset)[[group_var]] == group1)
        n_cells_group2 <- sum(colData(cell_subset)[[group_var]] == group2)
      }

      # Skip if not enough cells
      if (n_cells_group1 < min_cells || n_cells_group2 < min_cells) {
        message(paste0("Skipping ", cell_type, ": ", group1, " vs ", group2,
                       " (insufficient cells: ", n_cells_group1, ", ", n_cells_group2, ")"))
        next
      }

      # Run differential expression
      if (is_seurat) {
        Idents(cell_subset) <- cell_subset@meta.data[[group_var]]

        # [REVISED] Get all genes with no logFC threshold
        deg <- FindMarkers(cell_subset,
                           ident.1 = group1,
                           ident.2 = group2,
                           test.use = test.use,
                           logfc.threshold = 0,  # [REVISED] No threshold to keep all genes
                           min.pct = 0)  # [REVISED] Include all genes

      } else {
        # For SingleCellExperiment, use scran functions
        if (!requireNamespace("scran", quietly = TRUE)) {
          stop("Package 'scran' needed for SingleCellExperiment objects. Please install it.")
        }

        # Create factor for comparison
        group_factor <- factor(colData(cell_subset)[[group_var]])
        group_factor <- droplevels(group_factor)

        # [REVISED] Run findMarkers with no lfc threshold
        design <- model.matrix(~group_factor)
        deg <- scran::findMarkers(cell_subset,
                                  groups = group_factor,
                                  design = design,
                                  direction = "any",
                                  lfc = 0)  # [REVISED] No threshold

        # Format results to match Seurat output
        deg <- as.data.frame(deg[[group1]])
        deg$p_val_adj <- deg$FDR
        deg$avg_log2FC <- deg$logFC
      }

      # Add gene names as column
      deg$gene <- rownames(deg)

      # [REVISED] Add threshold passing info - keeping all genes but flagging those that pass threshold
      deg$passes_lfc_threshold <- abs(deg$avg_log2FC) >= logfc.threshold

      # Apply significance threshold - this defines "significant" based on p-value only
      deg$significant <- deg$p_val_adj < p_val_adj_threshold

      # Classify as up/down regulated
      deg$direction <- ifelse(deg$avg_log2FC > 0, "up", "down")

      # Store results
      comparison_name <- paste(group1, "vs", group2, sep = "_")
      result_key <- paste(cell_type, comparison_name, sep = "__")
      deg_results[[result_key]] <- deg

      # [REVISED] Compute summary statistics - only count genes that pass both p-value AND logFC threshold as DEGs
      n_up <- sum(deg$significant & deg$direction == "up" & deg$passes_lfc_threshold)
      n_down <- sum(deg$significant & deg$direction == "down" & deg$passes_lfc_threshold)

      # Add to summary data frame
      summary_row <- data.frame(
        cell_type = cell_type,
        comparison = comparison_name,
        group1 = group1,
        group2 = group2,
        n_cells_group1 = n_cells_group1,
        n_cells_group2 = n_cells_group2,
        n_genes_tested = nrow(deg),
        n_significant = sum(deg$significant & deg$passes_lfc_threshold), # [REVISED]
        n_up = n_up,
        n_down = n_down
      )

      summary_stats <- rbind(summary_stats, summary_row)
    }
  }

  # Return results
  results <- list(
    degs = deg_results,
    summary = summary_stats,
    metadata = list(
      cell_types = unique_cell_types,
      groups = unique_groups,
      comparisons = all_comparisons,
      parameters = list(
        min_cells = min_cells,
        test.use = test.use,
        logfc.threshold = logfc.threshold,
        p_val_adj_threshold = p_val_adj_threshold
      )
    )
  )

  return(results)
}
