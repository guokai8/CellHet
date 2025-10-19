#' Calculate Pseudotime Using Monocle3
#'
#' Wrapper function to calculate pseudotime trajectories using Monocle3.
#' Converts Seurat or SingleCellExperiment objects to CDS format, learns trajectory,
#' and returns pseudotime values.
#'
#' @param object Seurat or SingleCellExperiment object
#' @param reduction Dimensionality reduction to use (default: "pca" for SCE, "umap" for Seurat)
#' @param cluster_method Clustering method: "louvain" or "leiden" (default: "louvain")
#' @param root_cells Optional vector of cell names to use as trajectory root
#' @param root_type Optional cell type to use as trajectory root (requires cell_type_var)
#' @param cell_type_var Name of metadata column containing cell type annotations
#' @param num_dim Number of dimensions to use (default: 30)
#' @param partition_cells Whether to partition cells before learning graph (default: TRUE)
#' @param learn_graph_control Optional list of parameters for learn_graph
#' @param return_cds Whether to return the full CDS object (default: FALSE)
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return By default, returns the input object with pseudotime added to metadata.
#'   If return_cds=TRUE, returns a list with both the original object and CDS object.
#'
#' @details
#' This function provides a simplified interface to Monocle3 for pseudotime calculation.
#' It handles the conversion between object types, runs the trajectory analysis,
#' and adds the pseudotime values back to the original object's metadata.
#'
#' The pseudotime calculation involves:
#' 1. Converting to cell_data_set (CDS) format
#' 2. Preprocessing and reducing dimensions
#' 3. Clustering cells
#' 4. Learning the trajectory graph
#' 5. Ordering cells along the trajectory
#'
#' @examples
#' \dontrun{
#' # Basic usage with Seurat object
#' seurat_obj <- calculatePseudotime(seurat_obj)
#'
#' # Specify root cell type
#' seurat_obj <- calculatePseudotime(
#'   seurat_obj,
#'   root_type = "Stem_Cells",
#'   cell_type_var = "cell_type"
#' )
#'
#' # Return both original and CDS objects
#' results <- calculatePseudotime(seurat_obj, return_cds = TRUE)
#' seurat_obj <- results$object
#' cds <- results$cds
#'
#' # Use with SingleCellExperiment
#' sce <- calculatePseudotime(sce, reduction = "PCA")
#' }
#'
#' @export
calculatePseudotime <- function(object,
                               reduction = NULL,
                               cluster_method = "louvain",
                               root_cells = NULL,
                               root_type = NULL,
                               cell_type_var = NULL,
                               num_dim = 30,
                               partition_cells = TRUE,
                               learn_graph_control = NULL,
                               return_cds = FALSE,
                               verbose = TRUE) {

  # Check for monocle3
  if (!requireNamespace("monocle3", quietly = TRUE)) {
    stop("Package 'monocle3' is required for pseudotime calculation.\n",
         "Install with: BiocManager::install('monocle3')")
  }

  # Check object type
  is_seurat <- inherits(object, "Seurat")
  is_sce <- inherits(object, "SingleCellExperiment")

  if (!is_seurat && !is_sce) {
    stop("Input must be a Seurat or SingleCellExperiment object")
  }

  if (verbose) message("Starting pseudotime calculation with Monocle3...")

  # Set default reduction if not specified
  if (is.null(reduction)) {
    reduction <- if (is_seurat) "umap" else "PCA"
  }

  # Convert to CDS format
  if (verbose) message("Converting to CDS format...")

  if (is_seurat) {
    cds <- convertSeuratToCDS(object, reduction = reduction, verbose = verbose)
  } else {
    cds <- convertSCEToCDS(object, reduction = reduction, verbose = verbose)
  }

  # Cluster cells
  if (verbose) message("Clustering cells...")
  cds <- monocle3::cluster_cells(cds,
                                 reduction_method = reduction,
                                 cluster_method = cluster_method,
                                 num_dim = num_dim)

  # Partition cells
  if (partition_cells) {
    if (verbose) message("Partitioning cells...")
    cds <- monocle3::partition_cells(cds)
  }

  # Learn trajectory graph
  if (verbose) message("Learning trajectory graph...")
  if (is.null(learn_graph_control)) {
    cds <- monocle3::learn_graph(cds)
  } else {
    cds <- do.call(monocle3::learn_graph, c(list(cds), learn_graph_control))
  }

  # Determine root cells
  if (!is.null(root_type) && !is.null(cell_type_var)) {
    if (verbose) message(paste0("Using '", root_type, "' as trajectory root..."))

    # Get cells of root type
    if (is_seurat) {
      metadata <- object@meta.data
    } else {
      metadata <- as.data.frame(SummarizedExperiment::colData(object))
    }

    if (!cell_type_var %in% colnames(metadata)) {
      stop(paste0("Column '", cell_type_var, "' not found in metadata"))
    }

    root_cells <- rownames(metadata)[metadata[[cell_type_var]] == root_type]

    if (length(root_cells) == 0) {
      stop(paste0("No cells found with ", cell_type_var, " = '", root_type, "'"))
    }

    if (verbose) message(paste0("Found ", length(root_cells), " root cells"))
  }

  # Order cells (calculate pseudotime)
  if (verbose) message("Calculating pseudotime...")

  if (!is.null(root_cells)) {
    # Order with specified root
    cds <- monocle3::order_cells(cds, root_cells = root_cells)
  } else {
    # Interactive root selection would be needed here
    # For non-interactive use, we'll use the default
    warning("No root cells specified. Using default root selection.\n",
            "For better results, specify root_cells or root_type.")
    cds <- monocle3::order_cells(cds)
  }

  # Extract pseudotime values
  pseudotime_values <- monocle3::pseudotime(cds)

  # Add pseudotime to original object
  if (verbose) message("Adding pseudotime to object metadata...")

  if (is_seurat) {
    object$pseudotime <- pseudotime_values[colnames(object)]
    object$monocle3_partition <- cds@clusters$UMAP$partitions[colnames(object)]
    object$monocle3_cluster <- cds@clusters$UMAP$clusters[colnames(object)]
  } else {
    SummarizedExperiment::colData(object)$pseudotime <- pseudotime_values[colnames(object)]
    SummarizedExperiment::colData(object)$monocle3_partition <- cds@clusters[[1]]$partitions[colnames(object)]
    SummarizedExperiment::colData(object)$monocle3_cluster <- cds@clusters[[1]]$clusters[colnames(object)]
  }

  if (verbose) {
    message("Pseudotime calculation complete!")
    message(paste0("Added columns: pseudotime, monocle3_partition, monocle3_cluster"))
    message(paste0("Range: ", round(min(pseudotime_values, na.rm = TRUE), 2),
                  " to ", round(max(pseudotime_values, na.rm = TRUE), 2)))
  }

  # Return results
  if (return_cds) {
    return(list(
      object = object,
      cds = cds
    ))
  } else {
    return(object)
  }
}


#' Convert Seurat Object to Cell Data Set (CDS)
#'
#' @param seurat_obj Seurat object
#' @param reduction Reduction to use
#' @param verbose Whether to print messages
#'
#' @return A cell_data_set object
#' @keywords internal
convertSeuratToCDS <- function(seurat_obj, reduction = "umap", verbose = TRUE) {

  # Check if reduction exists
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(paste0("Reduction '", reduction, "' not found in Seurat object.\n",
               "Available reductions: ", paste(names(seurat_obj@reductions), collapse = ", ")))
  }

  # Get expression data
  if ("RNA" %in% names(seurat_obj@assays)) {
    expression_matrix <- Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  } else {
    expression_matrix <- Seurat::GetAssayData(seurat_obj, slot = "counts")
  }

  # Get cell metadata
  cell_metadata <- seurat_obj@meta.data

  # Get gene metadata
  gene_metadata <- data.frame(
    gene_short_name = rownames(expression_matrix),
    row.names = rownames(expression_matrix)
  )

  # Create CDS object
  cds <- monocle3::new_cell_data_set(
    expression_data = expression_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )

  # Add reduction
  SingleCellExperiment::reducedDims(cds)[[toupper(reduction)]] <- seurat_obj@reductions[[reduction]]@cell.embeddings

  # Add size factors if available
  if ("nCount_RNA" %in% colnames(cell_metadata)) {
    size_factors_vec <- cell_metadata$nCount_RNA / mean(cell_metadata$nCount_RNA)
    names(size_factors_vec) <- colnames(cds)
    BiocGenerics::sizeFactors(cds) <- size_factors_vec
  }

  return(cds)
}


#' Convert SingleCellExperiment to Cell Data Set (CDS)
#'
#' @param sce SingleCellExperiment object
#' @param reduction Reduction to use
#' @param verbose Whether to print messages
#'
#' @return A cell_data_set object
#' @keywords internal
convertSCEToCDS <- function(sce, reduction = "PCA", verbose = TRUE) {

  # Check if reduction exists
  if (!reduction %in% SingleCellExperiment::reducedDimNames(sce)) {
    stop(paste0("Reduction '", reduction, "' not found in SingleCellExperiment object.\n",
               "Available reductions: ",
               paste(SingleCellExperiment::reducedDimNames(sce), collapse = ", ")))
  }

  # Get expression data
  expression_matrix <- SummarizedExperiment::assay(sce, "counts")

  # Get cell metadata
  cell_metadata <- as.data.frame(SummarizedExperiment::colData(sce))

  # Get gene metadata
  gene_metadata <- as.data.frame(SummarizedExperiment::rowData(sce))
  if (nrow(gene_metadata) == 0) {
    gene_metadata <- data.frame(
      gene_short_name = rownames(expression_matrix),
      row.names = rownames(expression_matrix)
    )
  }

  # Create CDS object
  cds <- monocle3::new_cell_data_set(
    expression_data = expression_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )

  # Add reduction
  SingleCellExperiment::reducedDims(cds)[[toupper(reduction)]] <- SingleCellExperiment::reducedDim(sce, reduction)

  return(cds)
}


#' Plot Pseudotime Trajectory
#'
#' Create visualizations of pseudotime trajectories
#'
#' @param object Seurat or SingleCellExperiment object with pseudotime calculated
#' @param color_by Variable to color cells by (default: "pseudotime")
#' @param reduction Dimensionality reduction for plotting (default: "umap")
#' @param point_size Size of points (default: 0.5)
#' @param label_groups Whether to label groups (default: FALSE)
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' # After calculating pseudotime
#' plotPseudotime(seurat_obj)
#' plotPseudotime(seurat_obj, color_by = "cell_type")
#' }
#'
#' @export
plotPseudotime <- function(object,
                          color_by = "pseudotime",
                          reduction = "umap",
                          point_size = 0.5,
                          label_groups = FALSE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function to work. Please install it.")
  }

  # Check object type
  is_seurat <- inherits(object, "Seurat")
  is_sce <- inherits(object, "SingleCellExperiment")

  if (!is_seurat && !is_sce) {
    stop("Input must be a Seurat or SingleCellExperiment object")
  }

  # Extract data for plotting
  if (is_seurat) {
    # Check if pseudotime exists
    if (!"pseudotime" %in% colnames(object@meta.data)) {
      stop("Pseudotime not found. Run calculatePseudotime() first.")
    }

    # Get reduction coordinates
    if (!reduction %in% names(object@reductions)) {
      stop(paste0("Reduction '", reduction, "' not found"))
    }

    coords <- as.data.frame(object@reductions[[reduction]]@cell.embeddings[, 1:2])
    colnames(coords) <- c("Dim1", "Dim2")

    # Get metadata
    plot_data <- cbind(coords, object@meta.data)

  } else {
    # SCE object
    if (!"pseudotime" %in% colnames(SummarizedExperiment::colData(object))) {
      stop("Pseudotime not found. Run calculatePseudotime() first.")
    }

    # Get reduction coordinates
    if (!reduction %in% SingleCellExperiment::reducedDimNames(object)) {
      stop(paste0("Reduction '", reduction, "' not found"))
    }

    coords <- as.data.frame(SingleCellExperiment::reducedDim(object, reduction)[, 1:2])
    colnames(coords) <- c("Dim1", "Dim2")

    # Get metadata
    plot_data <- cbind(coords, as.data.frame(SummarizedExperiment::colData(object)))
  }

  # Check if color_by variable exists
  if (!color_by %in% colnames(plot_data)) {
    stop(paste0("Variable '", color_by, "' not found in metadata"))
  }

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Dim1, y = Dim2)) +
    ggplot2::geom_point(ggplot2::aes(color = .data[[color_by]]), size = point_size) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = paste0(toupper(reduction), "_1"),
      y = paste0(toupper(reduction), "_2"),
      title = "Pseudotime Trajectory",
      color = color_by
    ) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # Color scale based on variable type
  if (is.numeric(plot_data[[color_by]])) {
    p <- p + ggplot2::scale_color_viridis_c(option = "C")
  }

  return(p)
}


#' Get Pseudotime Statistics
#'
#' Calculate summary statistics for pseudotime values
#'
#' @param object Seurat or SingleCellExperiment object with pseudotime
#' @param group_by Optional variable to group statistics by
#'
#' @return A data frame with pseudotime statistics
#'
#' @examples
#' \dontrun{
#' # Overall statistics
#' getPseudotimeStats(seurat_obj)
#'
#' # Statistics by cell type
#' getPseudotimeStats(seurat_obj, group_by = "cell_type")
#' }
#'
#' @export
getPseudotimeStats <- function(object, group_by = NULL) {

  # Check object type
  is_seurat <- inherits(object, "Seurat")
  is_sce <- inherits(object, "SingleCellExperiment")

  if (!is_seurat && !is_sce) {
    stop("Input must be a Seurat or SingleCellExperiment object")
  }

  # Extract metadata
  if (is_seurat) {
    metadata <- object@meta.data
  } else {
    metadata <- as.data.frame(SummarizedExperiment::colData(object))
  }

  # Check if pseudotime exists
  if (!"pseudotime" %in% colnames(metadata)) {
    stop("Pseudotime not found. Run calculatePseudotime() first.")
  }

  # Calculate statistics
  if (is.null(group_by)) {
    # Overall statistics
    stats <- data.frame(
      n_cells = sum(!is.na(metadata$pseudotime)),
      n_na = sum(is.na(metadata$pseudotime)),
      mean = mean(metadata$pseudotime, na.rm = TRUE),
      median = median(metadata$pseudotime, na.rm = TRUE),
      min = min(metadata$pseudotime, na.rm = TRUE),
      max = max(metadata$pseudotime, na.rm = TRUE),
      sd = sd(metadata$pseudotime, na.rm = TRUE)
    )
  } else {
    # Grouped statistics
    if (!group_by %in% colnames(metadata)) {
      stop(paste0("Variable '", group_by, "' not found in metadata"))
    }

    stats <- do.call(rbind, lapply(unique(metadata[[group_by]]), function(grp) {
      grp_data <- metadata[metadata[[group_by]] == grp, ]
      data.frame(
        group = grp,
        n_cells = sum(!is.na(grp_data$pseudotime)),
        n_na = sum(is.na(grp_data$pseudotime)),
        mean = mean(grp_data$pseudotime, na.rm = TRUE),
        median = median(grp_data$pseudotime, na.rm = TRUE),
        min = min(grp_data$pseudotime, na.rm = TRUE),
        max = max(grp_data$pseudotime, na.rm = TRUE),
        sd = sd(grp_data$pseudotime, na.rm = TRUE)
      )
    }))
    colnames(stats)[1] <- group_by
  }

  return(stats)
}
