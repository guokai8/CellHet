#' Convert DEG Results from Seurat's FindMarkers to CellHet Format
#'
#' This function converts differential gene expression (DEG) results from Seurat's FindMarkers
#' format to the format expected by CellHet functions for visualization and analysis.
#'
#' @param deg_list A list of data frames, where each data frame contains DEG results from FindMarkers
#' @param p_val_threshold P-value threshold for significance (default: 0.05)
#' @param p_val_adj_threshold Adjusted p-value threshold for significance (default: 0.05)
#' @param logfc_threshold Log fold-change threshold for DEG identification (default: 0.25)
#' @param gene_col Column name for gene identifiers if they're not rownames (default: NULL)
#' @param add_rownames_as_genes Whether to add rownames as gene column if gene_col is NULL (default: TRUE)
#'
#' @return A list containing DEG results in CellHet format
#' @export
convertToCellHetDEG <- function(deg_list,
                                p_val_threshold = 0.05,
                                p_val_adj_threshold = 0.05,
                                logfc_threshold = 0.25,
                                gene_col = NULL,
                                add_rownames_as_genes = TRUE) {

  # Validate inputs
  if (!is.list(deg_list)) {
    stop("Input must be a list of data frames")
  }

  # Check if gene_col is provided or if rownames should be used
  if (is.null(gene_col) && !add_rownames_as_genes) {
    stop("Either gene_col must be provided or add_rownames_as_genes must be TRUE")
  }

  # Initialize CellHet style result structure
  cellhet_results <- list(
    degs = list(),
    summary = data.frame()
  )

  # Process each comparison
  for (comparison_name in names(deg_list)) {
    # Extract DEG data frame
    deg_data <- deg_list[[comparison_name]]

    # Skip if empty
    if (nrow(deg_data) == 0) {
      next
    }

    # Check if essential columns exist
    required_cols <- c("p_val_adj", "avg_log2FC")
    if (!all(required_cols %in% colnames(deg_data))) {
      warning(paste0("Skipping '", comparison_name, "': Missing required columns"))
      next
    }

    # Extract cell type and comparison from comparison_name
    # Format is often: CellType_GroupA_vs_GroupB
    # Handle the specific case in your data first
    if (grepl("_KO_vs_WT$", comparison_name)) {
      # In your data, format appears to be "CellType_KO_vs_WT"
      cell_type <- sub("_KO_vs_WT$", "", comparison_name)
      group1 <- "KO"
      group2 <- "WT"
      standard_comparison <- "KO_vs_WT"
    } else if (grepl("_WT_vs_KO$", comparison_name)) {
      # Handle reverse case too
      cell_type <- sub("_WT_vs_KO$", "", comparison_name)
      group1 <- "WT"
      group2 <- "KO"
      standard_comparison <- "WT_vs_KO"
    } else if (grepl("_vs_", comparison_name)) {
      # Generic handling for other "_vs_" patterns
      parts <- strsplit(comparison_name, "_vs_")[[1]]
      if (length(parts) >= 2) {
        # Last part is group2
        group2 <- parts[length(parts)]

        # Handle the first part which may contain cell type and group1
        first_part <- parts[1]
        # Try to extract group1 from the end of first_part if it ends with _KO, _WT, etc.
        if (grepl("_[A-Z]+$", first_part)) {
          # Extract the suffix that might be a group identifier
          potential_group <- sub(".*_([A-Z]+)$", "\\1", first_part)
          if (nchar(potential_group) <= 3) {  # Typical group IDs are short
            group1 <- potential_group
            cell_type <- sub(paste0("_", group1, "$"), "", first_part)
          } else {
            # If no clear group identifier, treat the whole first part as cell type
            cell_type <- first_part
            group1 <- "Group1"  # Default
          }
        } else {
          cell_type <- first_part
          group1 <- "Group1"  # Default
        }

        standard_comparison <- paste(group1, "vs", group2, sep="_")
      } else {
        # Fallback if pattern doesn't match
        cell_type <- comparison_name
        group1 <- "Group1"
        group2 <- "Group2"
        standard_comparison <- "Group1_vs_Group2"
      }
    } else {
      # Last resort fallback
      cell_type <- comparison_name
      group1 <- "Group1"
      group2 <- "Group2"
      standard_comparison <- "Group1_vs_Group2"
    }

    # Add gene column if needed
    if (!is.null(gene_col)) {
      # Use provided gene column
      if (!gene_col %in% colnames(deg_data)) {
        warning(paste0("Gene column '", gene_col, "' not found in '", comparison_name, "'. Using rownames."))
        deg_data$gene <- rownames(deg_data)
      } else {
        # Rename to 'gene' if it exists but has a different name
        if (gene_col != "gene") {
          deg_data$gene <- deg_data[[gene_col]]
        }
      }
    } else if (add_rownames_as_genes) {
      # Add rownames as gene column
      deg_data$gene <- rownames(deg_data)
    }

    # Add significance based on adjusted p-value
    deg_data$significant <- deg_data$p_val_adj < p_val_adj_threshold

    # Add logFC threshold flag
    deg_data$passes_lfc_threshold <- abs(deg_data$avg_log2FC) >= logfc_threshold

    # Add direction
    deg_data$direction <- ifelse(deg_data$avg_log2FC > 0, "up", "down")

    # Create cell counts if they don't exist
    if (!("n_cells_group1" %in% colnames(deg_data))) {
      deg_data$n_cells_group1 <- NA
    }
    if (!("n_cells_group2" %in% colnames(deg_data))) {
      deg_data$n_cells_group2 <- NA
    }

    # Create result key
    result_key <- paste(cell_type, standard_comparison, sep="__")

    # Store results
    cellhet_results$degs[[result_key]] <- deg_data

    # Calculate summary statistics
    n_up <- sum(deg_data$significant & deg_data$direction == "up" & deg_data$passes_lfc_threshold)
    n_down <- sum(deg_data$significant & deg_data$direction == "down" & deg_data$passes_lfc_threshold)
    n_significant <- sum(deg_data$significant & deg_data$passes_lfc_threshold)

    # Add to summary
    summary_row <- data.frame(
      cell_type = cell_type,
      comparison = standard_comparison,
      group1 = group1,
      group2 = group2,
      n_cells_group1 = NA,  # Will be filled if available
      n_cells_group2 = NA,  # Will be filled if available
      n_genes_tested = nrow(deg_data),
      n_significant = n_significant,
      n_up = n_up,
      n_down = n_down
    )

    # Check if we have cell counts
    if ("n_cells_group1" %in% colnames(deg_data) && !all(is.na(deg_data$n_cells_group1))) {
      summary_row$n_cells_group1 <- deg_data$n_cells_group1[1]
    }
    if ("n_cells_group2" %in% colnames(deg_data) && !all(is.na(deg_data$n_cells_group2))) {
      summary_row$n_cells_group2 <- deg_data$n_cells_group2[1]
    }

    # Add to summary data frame
    cellhet_results$summary <- rbind(cellhet_results$summary, summary_row)
  }

  # Add metadata to results
  cellhet_results$metadata <- list(
    cell_types = unique(cellhet_results$summary$cell_type),
    groups = unique(c(cellhet_results$summary$group1, cellhet_results$summary$group2)),
    comparisons = unique(cellhet_results$summary$comparison),
    parameters = list(
      min_cells = NA,
      test.use = NA,
      logfc.threshold = logfc_threshold,
      p_val_adj_threshold = p_val_threshold
    )
  )

  # Create all_degs data frame
  all_degs <- data.frame()

  for (result_key in names(cellhet_results$degs)) {
    # Extract cell type and comparison from key
    key_parts <- strsplit(result_key, "__")[[1]]
    cell_type <- key_parts[1]
    comparison <- key_parts[2]

    # Get DEG results
    degs <- cellhet_results$degs[[result_key]]

    # Add cell type and comparison columns if not present
    if (!("cell_type" %in% colnames(degs))) {
      degs$cell_type <- cell_type
    }
    if (!("comparison" %in% colnames(degs))) {
      degs$comparison <- comparison
    }

    # Split comparison into groups if not present
    if (!all(c("group1", "group2") %in% colnames(degs))) {
      comp_parts <- strsplit(comparison, "_vs_")[[1]]
      degs$group1 <- comp_parts[1]
      degs$group2 <- comp_parts[2]
    }

    # Append to combined data frame
    all_degs <- rbind(all_degs, degs)
  }

  # Add all_degs to results
  cellhet_results$all_degs <- all_degs

  return(cellhet_results)
}

