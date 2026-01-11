#' Compare DEGs Across Cell Types and Conditions
#'
#' @param object SingleCellExperiment or Seurat object with annotated cell types
#' @param group_var Name of the column in metadata containing group information
#' @param cell_type_var Name of the column in metadata containing cell type annotations
#' @param reference_group Optional reference group for comparisons (if NULL, all pairwise comparisons)
#' @param min_cells_per_group Minimum number of cells required per group for comparison
#' @param test_method Statistical test to use for differential expression
#' @param logfc_threshold Log fold-change threshold for DEGs identification
#' @param p_val_threshold P-value threshold for significance (unadjusted, default: NULL)
#' @param p_val_adj_threshold Adjusted p-value threshold for significance (default: 0.05)
#' @param custom_comparisons Optional list of custom comparisons to perform. Overrides reference_group if provided.
#'   Supports three formats: (1) 2-element vector: \code{c("GroupA", "GroupB")} for one-vs-one,
#'   (2) Named list: \code{list(test = c("Grp1", "Grp2"), reference = c("Grp3", "Grp4"))} for many-vs-many,
#'   (3) Vector with "vs": \code{c("Grp1", "Grp2", "vs", "Grp3", "Grp4")}. For multi-group comparisons,
#'   cells from the specified groups are pooled together for differential expression analysis.
#' @param cores Number of cores to use for parallel processing
#'
#' @return A list containing DEG results for each cell type and comparison
#'
#' @examples
#' \dontrun{
#' # Example 1: All pairwise comparisons (default)
#' results <- compareDEGs(seurat_obj, group_var = "condition", cell_type_var = "cell_type")
#'
#' # Example 2: Compare all groups against a reference
#' results <- compareDEGs(seurat_obj, reference_group = "Control")
#'
#' # Example 3: Custom one-vs-one comparisons
#' results <- compareDEGs(seurat_obj,
#'   custom_comparisons = list(c("Treatment1", "Control"), c("Treatment2", "Control")))
#'
#' # Example 4: Many-vs-many comparisons (pooling groups)
#' results <- compareDEGs(seurat_obj,
#'   custom_comparisons = list(
#'     list(test = c("Treatment1", "Treatment2"), reference = c("Control", "Baseline")),
#'     list(test = "Treatment3", reference = "Control")
#'   ))
#'
#' # Example 5: Using "vs" separator format
#' results <- compareDEGs(seurat_obj,
#'   custom_comparisons = list(c("TreatA", "TreatB", "vs", "Control", "Vehicle")))
#' }
#'
#' @export
compareDEGs <- function(object,
                        group_var = "condition",
                        cell_type_var = "cell_type",
                        reference_group = NULL,
                        min_cells_per_group = 3,
                        test_method = "wilcox",
                        logfc_threshold = 0.25,
                        p_val_threshold = NULL,
                        p_val_adj_threshold = 0.05,
                        custom_comparisons = NULL,
                        cores = 1) {

  # Check if input is Seurat or SingleCellExperiment
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

  # Extract unique cell types and groups
  unique_cell_types <- unique(metadata[[cell_type_var]])
  unique_groups <- unique(metadata[[group_var]])

  # Helper function to normalize comparison format
  normalize_comparison <- function(comp) {
    # Case 1: Named list with test/reference
    if (is.list(comp) && !is.null(names(comp)) && all(c("test", "reference") %in% names(comp))) {
      return(list(test = comp$test, reference = comp$reference))
    }

    # Case 2: Vector with "vs" separator
    if (is.character(comp) && "vs" %in% comp) {
      vs_idx <- which(comp == "vs")
      if (length(vs_idx) != 1) {
        stop("Comparison format error: only one 'vs' separator allowed")
      }
      if (vs_idx == 1 || vs_idx == length(comp)) {
        stop("Comparison format error: 'vs' separator cannot be at the beginning or end")
      }
      return(list(test = comp[1:(vs_idx-1)], reference = comp[(vs_idx+1):length(comp)]))
    }

    # Case 3: Legacy 2-element vector
    if (is.character(comp) && length(comp) == 2) {
      return(list(test = comp[1], reference = comp[2]))
    }

    # Case 4: Invalid format
    stop("Invalid comparison format. Must be: (1) 2-element vector c('A','B'), (2) named list list(test=c('A','B'), reference=c('C','D')), or (3) vector with 'vs' separator c('A','B','vs','C','D')")
  }

  # Helper function to create comparison name
  create_comparison_name <- function(test_groups, reference_groups) {
    test_str <- paste(test_groups, collapse = "+")
    ref_str <- paste(reference_groups, collapse = "+")
    return(paste(test_str, "vs", ref_str, sep = "_"))
  }

  # Define comparisons
  if (is.null(custom_comparisons)) {
    if (!is.null(reference_group)) {
      # Generate comparisons against reference group
      if (!reference_group %in% unique_groups) {
        stop(paste0("Reference group '", reference_group, "' not found in data"))
      }
      other_groups <- setdiff(unique_groups, reference_group)
      all_comparisons <- lapply(other_groups, function(g) list(test = g, reference = reference_group))
    } else {
      # Generate all pairwise comparisons
      pairs <- combn(unique_groups, 2, simplify = FALSE)
      all_comparisons <- lapply(pairs, function(p) list(test = p[1], reference = p[2]))
    }
  } else {
    # Normalize and validate user-defined comparisons
    all_comparisons <- lapply(custom_comparisons, normalize_comparison)

    # Validate all groups exist
    for (comp in all_comparisons) {
      invalid_test <- comp$test[!comp$test %in% unique_groups]
      invalid_ref <- comp$reference[!comp$reference %in% unique_groups]

      if (length(invalid_test) > 0) {
        stop(paste0("Test groups not found in data: ", paste(invalid_test, collapse = ", ")))
      }
      if (length(invalid_ref) > 0) {
        stop(paste0("Reference groups not found in data: ", paste(invalid_ref, collapse = ", ")))
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
      cell_subset <- subset(object, cells = which(metadata[[cell_type_var]] == cell_type))
    } else {
      cell_subset <- object[, metadata[[cell_type_var]] == cell_type]
    }

    # Process each comparison
    for (comp_idx in seq_along(all_comparisons)) {
      comparison <- all_comparisons[[comp_idx]]
      test_groups <- comparison$test
      reference_groups <- comparison$reference

      # Create comparison name
      comparison_name <- create_comparison_name(test_groups, reference_groups)

      # Check if enough cells for comparison (sum across all groups in each side)
      if (is_seurat) {
        n_cells_test <- sum(cell_subset@meta.data[[group_var]] %in% test_groups)
        n_cells_reference <- sum(cell_subset@meta.data[[group_var]] %in% reference_groups)
      } else {
        n_cells_test <- sum(colData(cell_subset)[[group_var]] %in% test_groups)
        n_cells_reference <- sum(colData(cell_subset)[[group_var]] %in% reference_groups)
      }

      # Skip if not enough cells
      if (n_cells_test < min_cells_per_group || n_cells_reference < min_cells_per_group) {
        message(paste0("Skipping ", cell_type, ": ", comparison_name,
                       " (insufficient cells: test=", n_cells_test, ", ref=", n_cells_reference, ")"))
        next
      }

      # Run differential expression
      if (is_seurat) {
        # Set identities using Seurat namespace
        Seurat::Idents(cell_subset) <- cell_subset@meta.data[[group_var]]

        # FindMarkers accepts vectors for both ident.1 and ident.2
        deg <- Seurat::FindMarkers(cell_subset,
                                   ident.1 = test_groups,      # Can be vector of identities
                                   ident.2 = reference_groups,  # Can be vector of identities
                                   test.use = test_method,
                                   logfc.threshold = 0,  # No threshold to keep all genes
                                   min.pct = 0)  # Include all genes

      } else {
        # For SingleCellExperiment, use scran functions
        if (!requireNamespace("scran", quietly = TRUE)) {
          stop("Package 'scran' needed for SingleCellExperiment objects. Please install it.")
        }

        # Create pooled group labels for many-vs-many comparisons
        cell_groups <- colData(cell_subset)[[group_var]]
        pooled_groups <- ifelse(
          cell_groups %in% test_groups,
          "TEST_POOL",
          ifelse(cell_groups %in% reference_groups, "REF_POOL", "OTHER")
        )

        # Filter to only test and reference cells
        keep_cells <- pooled_groups %in% c("TEST_POOL", "REF_POOL")
        cell_subset_filtered <- cell_subset[, keep_cells]
        pooled_groups <- factor(pooled_groups[keep_cells])

        # Run findMarkers
        design <- model.matrix(~pooled_groups)
        deg <- scran::findMarkers(cell_subset_filtered,
                                  groups = pooled_groups,
                                  design = design,
                                  direction = "any",
                                  lfc = 0)

        # Format results to match Seurat output
        deg <- as.data.frame(deg[["TEST_POOL"]])
        deg$p_val_adj <- deg$FDR
        deg$avg_log2FC <- deg$logFC
      }

      # Add gene names as column
      deg$gene <- rownames(deg)

      # Add threshold passing info - keeping all genes but flagging those that pass threshold
      deg$passes_lfc_threshold <- abs(deg$avg_log2FC) >= logfc_threshold

      # Apply significance threshold - use either p-value or adjusted p-value
      if (!is.null(p_val_threshold)) {
        # Use unadjusted p-value if specified
        deg$significant <- deg$p_val < p_val_threshold
      } else {
        # Use adjusted p-value (default)
        deg$significant <- deg$p_val_adj < p_val_adj_threshold
      }

      # Classify as up/down regulated
      deg$direction <- ifelse(deg$avg_log2FC > 0, "up", "down")

      # Add group information to each DEG result for downstream compatibility
      # For backward compatibility, use first group name for group1/group2
      deg$group1 <- test_groups[1]
      deg$group2 <- reference_groups[1]
      deg$n_cells_group1 <- n_cells_test
      deg$n_cells_group2 <- n_cells_reference

      # Add new columns for many-vs-many support
      deg$test_groups <- paste(test_groups, collapse = "+")
      deg$reference_groups <- paste(reference_groups, collapse = "+")

      # Store results
      result_key <- paste(cell_type, comparison_name, sep = "__")
      deg_results[[result_key]] <- deg

      # Compute summary statistics - only count genes that pass both p-value AND logFC threshold as DEGs
      n_up <- sum(deg$significant & deg$direction == "up" & deg$passes_lfc_threshold)
      n_down <- sum(deg$significant & deg$direction == "down" & deg$passes_lfc_threshold)

      # Add to summary data frame
      summary_row <- data.frame(
        cell_type = cell_type,
        comparison = comparison_name,
        test_groups = paste(test_groups, collapse = "+"),
        reference_groups = paste(reference_groups, collapse = "+"),
        n_cells_test = n_cells_test,
        n_cells_reference = n_cells_reference,
        n_genes_tested = nrow(deg),
        n_significant = sum(deg$significant & deg$passes_lfc_threshold),
        n_up = n_up,
        n_down = n_down,
        stringsAsFactors = FALSE
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
      reference_group = reference_group,
      parameters = list(
        min_cells_per_group = min_cells_per_group,
        test_method = test_method,
        logfc_threshold = logfc_threshold,
        p_val_threshold = p_val_threshold,
        p_val_adj_threshold = p_val_adj_threshold,
        using_adjusted_pval = is.null(p_val_threshold)
      )
    )
  )

  return(results)
}
