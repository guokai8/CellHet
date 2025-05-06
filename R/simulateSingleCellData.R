#' Generate Simulated Single-Cell Data
#'
#' Creates a simulated single-cell dataset with multiple cell types and conditions
#' for testing and demonstration purposes. Returns either a Seurat or SingleCellExperiment
#' object.
#'
#' @param n_cells Number of cells to simulate (default: 1000)
#' @param n_genes Number of genes to simulate (default: 2000)
#' @param n_cell_types Number of cell types to simulate (default: 5)
#' @param n_conditions Number of conditions to simulate (default: 2)
#' @param return_type Type of return object: "seurat" or "sce" (default: "seurat")
#' @param pct_de_genes Percentage of differentially expressed genes (default: 0.1)
#' @param effect_size Effect size for DE genes (default: 1.5)
#' @param batch_effect Whether to include batch effects (default: FALSE)
#' @param n_batches Number of batches if batch_effect is TRUE (default: 2)
#' @param batch_strength Strength of batch effect (default: 0.5)
#' @param add_trajectory Whether to add pseudotime trajectory (default: FALSE)
#' @param cell_type_names Optional vector of cell type names
#' @param condition_names Optional vector of condition names
#' @param seed Random seed for reproducibility (default: 42)
#'
#' @return A Seurat or SingleCellExperiment object with simulated data
#' @export
simulateSingleCellData <- function(n_cells = 1000,
                                   n_genes = 2000,
                                   n_cell_types = 5,
                                   n_conditions = 2,
                                   return_type = "seurat",
                                   pct_de_genes = 0.1,
                                   effect_size = 1.5,
                                   batch_effect = FALSE,
                                   n_batches = 2,
                                   batch_strength = 0.5,
                                   add_trajectory = FALSE,
                                   cell_type_names = NULL,
                                   condition_names = NULL,
                                   seed = 42) {

  # Check dependencies
  if (return_type == "seurat" && !requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' needed to return Seurat object. Please install it.")
  }
  if (return_type == "sce" && !requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' needed to return SCE object. Please install it.")
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Generate cell IDs - avoid using underscores as Seurat replaces them with dashes
  cell_ids <- paste0("cell", 1:n_cells)

  # Create metadata
  metadata <- data.frame(
    cell_id = cell_ids,
    row.names = cell_ids  # Important: set rownames to match cell IDs
  )

  # Assign cell types
  if (is.null(cell_type_names)) {
    cell_type_names <- paste0("CellType", 1:n_cell_types)
  } else if (length(cell_type_names) < n_cell_types) {
    warning("Not enough cell type names provided. Using default naming for remaining types.")
    cell_type_names <- c(cell_type_names, paste0("CellType", (length(cell_type_names)+1):n_cell_types))
  }

  cell_per_type <- n_cells / n_cell_types
  metadata$cell_type <- rep(cell_type_names, each = cell_per_type)

  # Assign conditions
  if (is.null(condition_names)) {
    condition_names <- paste0("Condition", 1:n_conditions)
  } else if (length(condition_names) < n_conditions) {
    warning("Not enough condition names provided. Using default naming for remaining conditions.")
    condition_names <- c(condition_names, paste0("Condition", (length(condition_names)+1):n_conditions))
  }

  cells_per_type_condition <- cell_per_type / n_conditions
  metadata$condition <- rep(rep(condition_names, each = cells_per_type_condition), n_cell_types)

  # Add batch information if requested
  if (batch_effect) {
    metadata$batch <- sample(paste0("Batch", 1:n_batches), n_cells, replace = TRUE)
  }

  # Generate gene names - avoid underscores
  gene_names <- paste0("gene", 1:n_genes)

  # Create count matrix - ensure colnames match rownames of metadata
  counts <- matrix(0, nrow = n_genes, ncol = n_cells)
  rownames(counts) <- gene_names
  colnames(counts) <- cell_ids  # Match cell IDs

  # Generate base expression profiles for each cell type
  cell_type_profiles <- matrix(0, nrow = n_genes, ncol = n_cell_types)

  # Set base expression - some genes specific to each cell type
  genes_per_cell_type <- round(n_genes * 0.05)  # 5% of genes specific to each cell type

  for (i in 1:n_cell_types) {
    # Base expression for all genes
    cell_type_profiles[, i] <- rnbinom(n_genes, size = 1, prob = 0.7)

    # Marker genes for this cell type (higher expression)
    marker_genes <- ((i-1) * genes_per_cell_type + 1):(i * genes_per_cell_type)
    cell_type_profiles[marker_genes, i] <- cell_type_profiles[marker_genes, i] * 5
  }

  # Fill count matrix based on cell type profiles
  for (i in 1:n_cells) {
    cell_type_idx <- match(metadata$cell_type[i], cell_type_names)
    base_expression <- cell_type_profiles[, cell_type_idx]

    # Add noise
    counts[, i] <- rnbinom(n_genes, mu = base_expression, size = 10)
  }

  # Add differential expression between conditions
  n_de_genes <- round(n_genes * pct_de_genes)
  de_genes <- sample(gene_names, n_de_genes)

  # For each cell type, apply DE effects
  for (ct in cell_type_names) {
    for (cond in condition_names[-1]) {  # Skip first condition (reference)
      # Get cells of this type and condition
      cells_idx <- which(metadata$cell_type == ct & metadata$condition == cond)

      # Apply effect size to DE genes
      effect_direction <- sample(c(-1, 1), n_de_genes, replace = TRUE)  # Up or down regulation
      effect_multiplier <- 1 + (effect_size * effect_direction)

      for (i in 1:length(de_genes)) {
        gene_idx <- match(de_genes[i], gene_names)
        counts[gene_idx, cells_idx] <- round(counts[gene_idx, cells_idx] * effect_multiplier[i])
      }
    }
  }

  # Apply batch effects if requested
  if (batch_effect) {
    for (batch in paste0("Batch", 1:n_batches)) {
      # Get cells in this batch
      batch_cells <- which(metadata$batch == batch)

      # Skip first batch (reference)
      if (batch == "Batch1") next

      # Select random genes for batch effect
      batch_genes <- sample(gene_names, round(n_genes * 0.2))  # 20% of genes affected by batch

      # Apply batch effect
      for (gene in batch_genes) {
        gene_idx <- match(gene, gene_names)
        batch_multiplier <- 1 + rnorm(1, mean = 0, sd = batch_strength)
        counts[gene_idx, batch_cells] <- round(counts[gene_idx, batch_cells] * batch_multiplier)
      }
    }
  }

  # Add trajectory if requested
  if (add_trajectory) {
    # Generate pseudotime values
    metadata$pseudotime <- rep(0, n_cells)

    # Different trajectory patterns for different cell types
    for (ct_idx in 1:n_cell_types) {
      ct <- cell_type_names[ct_idx]
      ct_cells <- which(metadata$cell_type == ct)

      # Generate random trajectory values with some pattern
      metadata$pseudotime[ct_cells] <- sort(runif(length(ct_cells))) + (ct_idx - 1) * 0.2
    }

    # Add genes that change along trajectory
    trajectory_genes <- sample(setdiff(gene_names, de_genes), round(n_genes * 0.05))

    for (gene in trajectory_genes) {
      gene_idx <- match(gene, gene_names)

      # Apply different patterns to different genes
      pattern_type <- sample(1:3, 1)

      if (pattern_type == 1) {
        # Linear increase
        trajectory_effect <- metadata$pseudotime
      } else if (pattern_type == 2) {
        # Quadratic
        trajectory_effect <- metadata$pseudotime^2
      } else {
        # Sigmoid
        trajectory_effect <- 1 / (1 + exp(-5 * (metadata$pseudotime - 0.5)))
      }

      # Scale effect
      trajectory_effect <- trajectory_effect / max(trajectory_effect) * 3

      # Apply effect
      counts[gene_idx, ] <- counts[gene_idx, ] * (1 + trajectory_effect)
    }
  }

  # Ensure non-negative counts
  counts[counts < 0] <- 0

  # Round all counts to integers
  counts <- round(counts)

  # Create return object
  if (return_type == "seurat") {
    # Create Seurat object - make sure metadata rownames match counts colnames
    seurat_obj <- Seurat::CreateSeuratObject(counts = counts)

    # Now add our custom metadata
    for (col in colnames(metadata)) {
      if (col != "cell_id") {  # Skip cell_id as it's redundant with rownames
        seurat_obj[[col]] <- metadata[[col]]
      }
    }

    # Normalize and scale data
    seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
    seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, nfeatures = 1000, verbose = FALSE)
    seurat_obj <- Seurat::ScaleData(seurat_obj, verbose = FALSE)

    # Run PCA and UMAP
    seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
    seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)

    return(seurat_obj)

  } else if (return_type == "sce") {
    # Create SingleCellExperiment object
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("Package 'SingleCellExperiment' needed for SCE object. Please install it.")
    }

    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = counts),
      colData = metadata
    )

    # Add logcounts if requested
    if (requireNamespace("scater", quietly = TRUE)) {
      sce <- scater::logNormCounts(sce)

      # Add PCA and UMAP
      sce <- scater::runPCA(sce)
      sce <- scater::runUMAP(sce)
    }

    return(sce)

  } else {
    stop("return_type must be either 'seurat' or 'sce'")
  }
}
