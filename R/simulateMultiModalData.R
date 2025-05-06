#' Generate Simulated Multi-Modal Single-Cell Data
#'
#' Creates a simulated multi-modal single-cell dataset with multiple cell types
#' and conditions for testing integration functions. Returns either a list of
#' Seurat objects or SingleCellExperiment objects.
#'
#' @param n_cells Number of cells to simulate (default: 1000)
#' @param n_modalities Number of modalities to simulate (default: 2)
#' @param n_features Vector of features per modality (default: c(2000, 500))
#' @param n_cell_types Number of cell types to simulate (default: 5)
#' @param n_conditions Number of conditions to simulate (default: 2)
#' @param return_type Type of return object: "seurat" or "sce" (default: "seurat")
#' @param pct_de_features Percentage of differentially expressed features (default: 0.1)
#' @param effect_size Effect size for DE features (default: 1.5)
#' @param shared_cell_types Whether to use the same cell types across modalities (default: TRUE)
#' @param cell_overlap Percentage of cells shared between modalities (default: 0.7)
#' @param feature_overlap Percentage of features with correlated expression (default: 0.3)
#' @param correlation_strength Strength of correlation between shared features (default: 0.7)
#' @param batch_effect Whether to include batch effects (default: FALSE)
#' @param n_batches Number of batches if batch_effect is TRUE (default: 2)
#' @param modality_names Optional vector of modality names
#' @param cell_type_names Optional vector of cell type names
#' @param condition_names Optional vector of condition names
#' @param add_trajectory Whether to add pseudotime trajectory (default: FALSE)
#' @param seed Random seed for reproducibility (default: 42)
#'
#' @return A list of Seurat or SingleCellExperiment objects with simulated multi-modal data
#' @export
simulateMultiModalData <- function(n_cells = 1000,
                                   n_modalities = 2,
                                   n_features = c(2000, 500),
                                   n_cell_types = 5,
                                   n_conditions = 2,
                                   return_type = "seurat",
                                   pct_de_features = 0.1,
                                   effect_size = 1.5,
                                   shared_cell_types = TRUE,
                                   cell_overlap = 0.7,
                                   feature_overlap = 0.3,
                                   correlation_strength = 0.7,
                                   batch_effect = FALSE,
                                   n_batches = 2,
                                   modality_names = NULL,
                                   cell_type_names = NULL,
                                   condition_names = NULL,
                                   add_trajectory = FALSE,
                                   seed = 42) {

  # Check dependencies
  if (return_type == "seurat" && !requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' needed to return Seurat object. Please install it.")
  }
  if (return_type == "sce" && !requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' needed to return SCE object. Please install it.")
  }

  # Validate inputs
  if (length(n_features) != n_modalities) {
    message("Length of n_features doesn't match n_modalities. Using the first value for all modalities.")
    n_features <- rep(n_features[1], n_modalities)
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Set modality names
  if (is.null(modality_names)) {
    modality_names <- paste0("Modality", 1:n_modalities)
  } else if (length(modality_names) < n_modalities) {
    warning("Not enough modality names provided. Using default naming for remaining modalities.")
    modality_names <- c(modality_names, paste0("Modality", (length(modality_names)+1):n_modalities))
  }

  # Set cell type names
  if (is.null(cell_type_names)) {
    cell_type_names <- paste0("CellType", 1:n_cell_types)
  } else if (length(cell_type_names) < n_cell_types) {
    warning("Not enough cell type names provided. Using default naming for remaining types.")
    cell_type_names <- c(cell_type_names, paste0("CellType", (length(cell_type_names)+1):n_cell_types))
  }

  # Set condition names
  if (is.null(condition_names)) {
    condition_names <- paste0("Condition", 1:n_conditions)
  } else if (length(condition_names) < n_conditions) {
    warning("Not enough condition names provided. Using default naming for remaining conditions.")
    condition_names <- c(condition_names, paste0("Condition", (length(condition_names)+1):n_conditions))
  }

  # Generate cell metadata for each modality
  modality_data <- list()
  cell_ids_by_modality <- list()

  # Determine cell overlap
  total_cells <- round(n_cells / cell_overlap)
  # Use NO underscores in cell IDs
  cell_ids <- paste0("cell", 1:total_cells)

  for (m in 1:n_modalities) {
    # Select cells for this modality
    if (m == 1) {
      # First modality gets n_cells
      mod_cell_ids <- cell_ids[1:n_cells]
    } else {
      # Subsequent modalities get overlapping and unique cells
      n_overlap <- round(n_cells * cell_overlap)
      n_unique <- n_cells - n_overlap

      overlap_ids <- sample(cell_ids_by_modality[[1]], n_overlap)
      unique_ids <- setdiff(cell_ids, unlist(cell_ids_by_modality))[1:n_unique]

      mod_cell_ids <- c(overlap_ids, unique_ids)
    }

    # Store cell IDs for this modality
    cell_ids_by_modality[[m]] <- mod_cell_ids

    # Create cell metadata
    cells_per_type <- length(mod_cell_ids) / n_cell_types
    cells_per_type_condition <- cells_per_type / n_conditions

    metadata <- data.frame(
      cell_id = mod_cell_ids,
      row.names = mod_cell_ids  # Important: set rownames to match cell IDs
    )

    # Assign cell types - either same as modality 1 for shared cells or new assignments
    if (m == 1 || !shared_cell_types) {
      # Random assignment
      metadata$cell_type <- rep(cell_type_names, each = cells_per_type)

      # Shuffle cell types except for first modality
      if (m > 1) {
        metadata$cell_type <- sample(metadata$cell_type)
      }
    } else {
      # For shared cell types, maintain cell type assignment for overlapping cells
      metadata$cell_type <- NA

      # For overlapping cells, get cell type from modality 1
      overlap_cells <- intersect(mod_cell_ids, cell_ids_by_modality[[1]])
      overlap_indices_mod <- match(overlap_cells, mod_cell_ids)
      overlap_indices_mod1 <- match(overlap_cells, cell_ids_by_modality[[1]])

      metadata$cell_type[overlap_indices_mod] <- modality_data[[1]]$metadata$cell_type[overlap_indices_mod1]

      # For unique cells, assign random cell types
      unique_indices <- which(is.na(metadata$cell_type))
      metadata$cell_type[unique_indices] <- sample(cell_type_names, length(unique_indices), replace = TRUE)
    }

    # Assign conditions - can be different for each modality
    metadata$condition <- rep(rep(condition_names, each = cells_per_type_condition), n_cell_types)
    if (m > 1) {
      metadata$condition <- sample(metadata$condition)  # Shuffle for variety
    }

    # Add batch information if requested
    if (batch_effect) {
      metadata$batch <- sample(paste0("Batch", 1:n_batches), length(mod_cell_ids), replace = TRUE)
    }

    # Generate feature names - avoid underscores
    feature_names <- paste0(modality_names[m], "feature", 1:n_features[m])

    # Create count matrix - ensure colnames match rownames of metadata
    counts <- matrix(0, nrow = n_features[m], ncol = length(mod_cell_ids))
    rownames(counts) <- feature_names
    colnames(counts) <- mod_cell_ids

    # Generate base expression profiles for each cell type
    cell_type_profiles <- matrix(0, nrow = n_features[m], ncol = n_cell_types)

    # Set base expression - some features specific to each cell type
    features_per_cell_type <- round(n_features[m] * 0.05)  # 5% of features specific to each cell type

    for (i in 1:n_cell_types) {
      # Base expression for all features
      cell_type_profiles[, i] <- rnbinom(n_features[m], size = 1, prob = 0.7)

      # Marker features for this cell type (higher expression)
      marker_features <- ((i-1) * features_per_cell_type + 1):(i * features_per_cell_type)
      cell_type_profiles[marker_features, i] <- cell_type_profiles[marker_features, i] * 5
    }

    # Fill count matrix based on cell type profiles
    for (i in 1:length(mod_cell_ids)) {
      cell_type_idx <- match(metadata$cell_type[i], cell_type_names)
      base_expression <- cell_type_profiles[, cell_type_idx]

      # Add noise
      counts[, i] <- rnbinom(n_features[m], mu = base_expression, size = 10)
    }

    # Add differential expression between conditions
    n_de_features <- round(n_features[m] * pct_de_features)
    de_features <- sample(feature_names, n_de_features)

    # For each cell type, apply DE effects
    for (ct in cell_type_names) {
      for (cond in condition_names[-1]) {  # Skip first condition (reference)
        # Get cells of this type and condition
        cells_idx <- which(metadata$cell_type == ct & metadata$condition == cond)

        # Apply effect size to DE features
        effect_direction <- sample(c(-1, 1), n_de_features, replace = TRUE)  # Up or down regulation
        effect_multiplier <- 1 + (effect_size * effect_direction)

        for (i in 1:length(de_features)) {
          feature_idx <- match(de_features[i], feature_names)
          counts[feature_idx, cells_idx] <- round(counts[feature_idx, cells_idx] * effect_multiplier[i])
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

        # Select random features for batch effect
        batch_features <- sample(feature_names, round(n_features[m] * 0.2))  # 20% of features affected by batch

        # Apply batch effect
        for (feature in batch_features) {
          feature_idx <- match(feature, feature_names)
          batch_multiplier <- 1 + rnorm(1, mean = 0, sd = 0.5)
          counts[feature_idx, batch_cells] <- round(counts[feature_idx, batch_cells] * batch_multiplier)
        }
      }
    }

    # Add trajectory if requested
    if (add_trajectory) {
      # Generate pseudotime values - try to maintain consistent trajectory across modalities
      if (m == 1 || !any(grepl("pseudotime", colnames(modality_data[[1]]$metadata)))) {
        metadata$pseudotime <- rep(0, length(mod_cell_ids))

        # Different trajectory patterns for different cell types
        for (ct_idx in 1:n_cell_types) {
          ct <- cell_type_names[ct_idx]
          ct_cells <- which(metadata$cell_type == ct)

          # Generate random trajectory values with some pattern
          metadata$pseudotime[ct_cells] <- sort(runif(length(ct_cells))) + (ct_idx - 1) * 0.2
        }
      } else {
        # Use pseudotime from first modality for overlapping cells
        metadata$pseudotime <- rep(NA, length(mod_cell_ids))

        overlap_cells <- intersect(mod_cell_ids, cell_ids_by_modality[[1]])
        overlap_indices_mod <- match(overlap_cells, mod_cell_ids)
        overlap_indices_mod1 <- match(overlap_cells, cell_ids_by_modality[[1]])

        metadata$pseudotime[overlap_indices_mod] <- modality_data[[1]]$metadata$pseudotime[overlap_indices_mod1]

        # For new cells, generate values based on cell type
        for (ct in cell_type_names) {
          missing_cells <- which(is.na(metadata$pseudotime) & metadata$cell_type == ct)
          if (length(missing_cells) > 0) {
            # Get range from existing cells
            existing_range <- range(metadata$pseudotime[!is.na(metadata$pseudotime) &
                                                          metadata$cell_type == ct],
                                    na.rm = TRUE)

            if (all(is.na(existing_range)) || length(existing_range) == 0) {
              # No existing cells of this type, generate new values
              ct_idx <- match(ct, cell_type_names)
              metadata$pseudotime[missing_cells] <- runif(length(missing_cells)) + (ct_idx - 1) * 0.2
            } else {
              # Generate values in same range as existing cells
              metadata$pseudotime[missing_cells] <- runif(length(missing_cells),
                                                          min = existing_range[1],
                                                          max = existing_range[2])
            }
          }
        }
      }

      # Add features that change along trajectory
      trajectory_features <- sample(setdiff(feature_names, de_features), round(n_features[m] * 0.05))

      for (feature in trajectory_features) {
        feature_idx <- match(feature, feature_names)

        # Apply different patterns to different features
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
        counts[feature_idx, ] <- counts[feature_idx, ] * (1 + trajectory_effect)
      }
    }

    # Ensure non-negative counts
    counts[counts < 0] <- 0
    counts <- round(counts)

    # Store data for this modality
    modality_data[[m]] <- list(
      metadata = metadata,
      counts = counts,
      modality = modality_names[m]
    )
  }

  # Add correlated features between modalities
  if (n_modalities > 1 && feature_overlap > 0) {
    for (m in 2:n_modalities) {
      # Number of features to correlate
      n_corr_features <- round(min(n_features[1], n_features[m]) * feature_overlap)

      # Select features from both modalities
      mod1_features <- sample(rownames(modality_data[[1]]$counts), n_corr_features)
      mod_m_features <- sample(rownames(modality_data[[m]]$counts), n_corr_features)

      # Find overlapping cells
      overlap_cells <- intersect(colnames(modality_data[[1]]$counts),
                                 colnames(modality_data[[m]]$counts))

      # For each pair of features, add correlation
      for (i in 1:n_corr_features) {
        mod1_feature <- mod1_features[i]
        mod_m_feature <- mod_m_features[i]

        # Get expression in overlapping cells
        mod1_expr <- modality_data[[1]]$counts[mod1_feature, overlap_cells]

        # Create correlated expression
        noise <- rnorm(length(overlap_cells), mean = 0, sd = 1 - correlation_strength)
        mod_m_expr <- correlation_strength * mod1_expr + noise

        # Scale to match original distribution
        mod_m_expr <- mod_m_expr * (mean(modality_data[[m]]$counts[mod_m_feature, overlap_cells]) /
                                      mean(mod_m_expr))

        # Replace values
        modality_data[[m]]$counts[mod_m_feature, overlap_cells] <- round(pmax(0, mod_m_expr))
      }
    }
  }

  # Create return objects
  return_objects <- list()

  for (m in 1:n_modalities) {
    if (return_type == "seurat") {
      # Create Seurat object
      seurat_obj <- Seurat::CreateSeuratObject(counts = modality_data[[m]]$counts)

      # Add metadata columns
      for (col in colnames(modality_data[[m]]$metadata)) {
        if (col != "cell_id") {  # Skip cell_id as it's redundant with rownames
          seurat_obj[[col]] <- modality_data[[m]]$metadata[[col]]
        }
      }

      # Add modality information
      seurat_obj$modality <- modality_names[m]

      # Normalize and scale data
      seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
      seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, nfeatures = min(1000, n_features[m]), verbose = FALSE)
      seurat_obj <- Seurat::ScaleData(seurat_obj, verbose = FALSE)

      # Run PCA and UMAP
      seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = min(30, n_features[m] - 1), verbose = FALSE)
      seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:min(30, n_features[m] - 1), verbose = FALSE)

      return_objects[[modality_names[m]]] <- seurat_obj

    } else if (return_type == "sce") {
      # Create SingleCellExperiment object
      sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = modality_data[[m]]$counts),
        colData = modality_data[[m]]$metadata
      )

      # Add modality information
      SummarizedExperiment::colData(sce)$modality <- modality_names[m]

      # Add logcounts if requested
      if (requireNamespace("scater", quietly = TRUE)) {
        sce <- scater::logNormCounts(sce)

        # Add PCA and UMAP
        sce <- scater::runPCA(sce, ncomponents = min(30, n_features[m] - 1))
        sce <- scater::runUMAP(sce)
      }

      return_objects[[modality_names[m]]] <- sce
    }
  }

  # Add metadata about simulation
  attr(return_objects, "simulation") <- list(
    n_cells = n_cells,
    n_modalities = n_modalities,
    n_features = n_features,
    n_cell_types = n_cell_types,
    n_conditions = n_conditions,
    shared_cell_types = shared_cell_types,
    cell_overlap = cell_overlap,
    feature_overlap = feature_overlap,
    batch_effect = batch_effect,
    add_trajectory = add_trajectory,
    modality_names = modality_names,
    cell_type_names = cell_type_names,
    condition_names = condition_names,
    seed = seed
  )

  return(return_objects)
}
