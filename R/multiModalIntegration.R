# multiModalIntegration.R

#' Integrate Multiple Modalities for Cell Heterogeneity Analysis
#'
#' Analyzes cellular heterogeneity across multiple data modalities (RNA, protein, chromatin, etc.)
#' to identify shared and modality-specific patterns.
#'
#' @param data_list List of data objects (Seurat or SingleCellExperiment) for different modalities
#' @param modality_names Names of modalities corresponding to each item in data_list
#' @param cell_type_var Name of the variable containing cell type annotations
#' @param integration_method Method for integration: "canonical_correlation", "mnn", "harmony", or "seurat" (default: "canonical_correlation")
#' @param dim_reduction Post-integration dimension reduction: "umap", "tsne", or "none" (default: "umap")
#' @param n_dims Number of dimensions to use for integration (default: 30)
#' @param n_features Number of features to use from each modality (default: 2000)
#' @param scale_data Logical, whether to scale data before integration (default: TRUE)
#' @param seed Random seed for reproducibility
#'
#' @return A list containing integrated data and analysis results
#' @export
integrateMultiModalData <- function(data_list,
                                    modality_names = NULL,
                                    cell_type_var = "cell_type",
                                    integration_method = "canonical_correlation",
                                    dim_reduction = "umap",
                                    n_dims = 30,
                                    n_features = 2000,
                                    scale_data = TRUE,
                                    seed = 42) {

  # Validate input
  if (!is.list(data_list)) {
    stop("data_list must be a list of data objects")
  }

  n_modalities <- length(data_list)
  if (n_modalities < 2) {
    stop("At least two modalities are required for integration")
  }

  # Set modality names if not provided
  if (is.null(modality_names)) {
    modality_names <- paste0("Modality", 1:n_modalities)
  } else if (length(modality_names) != n_modalities) {
    stop("Length of modality_names must match length of data_list")
  }

  # Determine object types and validate
  object_types <- sapply(data_list, function(x) {
    if (inherits(x, "Seurat")) {
      return("Seurat")
    } else if (inherits(x, "SingleCellExperiment")) {
      return("SCE")
    } else {
      stop("Data objects must be either Seurat or SingleCellExperiment")
    }
  })

  # Set seed for reproducibility
  set.seed(seed)

  # Initialize results list
  results <- list(
    metadata = list(
      modality_names = modality_names,
      object_types = object_types,
      integration_method = integration_method,
      dim_reduction = dim_reduction,
      n_dims = n_dims,
      n_features = n_features
    ),
    integration_results = NULL,
    integration_embeddings = NULL,
    plots = list()
  )

  # Perform integration based on method
  if (integration_method == "canonical_correlation") {
    # Canonical Correlation Analysis (CCA)

    # For Seurat objects
    if (all(object_types == "Seurat")) {
      if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' required for CCA integration")
      }

      # Each modality should have variable features identified
      for (i in 1:n_modalities) {
        if (length(Seurat::VariableFeatures(data_list[[i]])) == 0) {
          message(paste0("Finding variable features for ", modality_names[i]))
          data_list[[i]] <- Seurat::FindVariableFeatures(data_list[[i]],
                                                         selection.method = "vst",
                                                         nfeatures = n_features)
        }
      }

      # Perform CCA integration
      message("Performing CCA integration...")
      integration_anchors <- Seurat::FindIntegrationAnchors(
        object.list = data_list,
        dims = 1:n_dims,
        anchor.features = n_features,
        scale = scale_data
      )

      # Create integrated object
      integrated_data <- Seurat::IntegrateData(
        anchorset = integration_anchors,
        dims = 1:n_dims
      )

      # Scale integrated data
      integrated_data <- Seurat::ScaleData(integrated_data)

      # Run dimension reduction
      integrated_data <- Seurat::RunPCA(integrated_data, npcs = n_dims)

      if (dim_reduction == "umap") {
        integrated_data <- Seurat::RunUMAP(integrated_data, dims = 1:n_dims)
      } else if (dim_reduction == "tsne") {
        integrated_data <- Seurat::RunTSNE(integrated_data, dims = 1:n_dims)
      }

      # Extract integration embeddings
      if (dim_reduction == "umap") {
        embeddings <- Seurat::Embeddings(integrated_data, reduction = "umap")
      } else if (dim_reduction == "tsne") {
        embeddings <- Seurat::Embeddings(integrated_data, reduction = "tsne")
      } else {
        embeddings <- Seurat::Embeddings(integrated_data, reduction = "pca")[, 1:2]
      }

      # Store results
      results$integration_results <- integrated_data
      results$integration_embeddings <- embeddings

      # Create plots if ggplot2 is available
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        # Extract metadata
        metadata <- integrated_data@meta.data

        # Create dimension reduction plot
        if (cell_type_var %in% colnames(metadata)) {
          plot_data <- data.frame(
            x = embeddings[, 1],
            y = embeddings[, 2],
            cell_type = metadata[[cell_type_var]]
          )

          dim_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = cell_type)) +
            ggplot2::geom_point(size = 1, alpha = 0.7) +
            ggplot2::theme_minimal() +
            ggplot2::labs(
              title = paste0("Integrated ", dim_reduction, " Plot"),
              x = paste0(dim_reduction, "1"),
              y = paste0(dim_reduction, "2"),
              color = "Cell Type"
            )

          results$plots$dim_plot <- dim_plot
        }
      }

    } else if (all(object_types == "SCE")) {
      # For SingleCellExperiment objects

      if (!requireNamespace("batchelor", quietly = TRUE)) {
        stop("Package 'batchelor' required for integration of SingleCellExperiment objects")
      }

      # Prepare list for batchelor
      sce_list <- data_list

      # Apply dimensionality reduction to each modality
      reduced_dims <- lapply(sce_list, function(sce) {
        if (!requireNamespace("scater", quietly = TRUE)) {
          stop("Package 'scater' required for processing SingleCellExperiment objects")
        }
        # Run PCA
        sce <- scater::runPCA(sce, ncomponents = n_dims)
        return(SingleCellExperiment::reducedDim(sce, "PCA"))
      })

      # Perform multiCCA using batchelor
      message("Performing multiCCA integration...")
      integration_results <- batchelor::multiBatchPCA(
        sce_list,
        d = n_dims,
        weights = TRUE
      )

      # Correct for batch effects
      corrected <- batchelor::reducedMNN(
        integration_results$corrected,
        batch = rep(1:n_modalities, sapply(sce_list, ncol)),
        k = 20,
        d = n_dims
      )

      # Store corrected dimensions in a new SingleCellExperiment
      integrated_data <- SingleCellExperiment::SingleCellExperiment(
        assays = list(logcounts = matrix(0, nrow = 0, ncol = ncol(corrected$corrected)))
      )

      SingleCellExperiment::reducedDim(integrated_data, "MNN") <- corrected$corrected

      # Run dimension reduction on corrected data
      if (dim_reduction == "umap") {
        if (!requireNamespace("uwot", quietly = TRUE)) {
          warning("Package 'uwot' required for UMAP. Using PCA instead.")
          embeddings <- corrected$corrected[, 1:2]
        } else {
          umap_res <- uwot::umap(corrected$corrected, n_components = 2)
          embeddings <- umap_res
          SingleCellExperiment::reducedDim(integrated_data, "UMAP") <- umap_res
        }
      } else if (dim_reduction == "tsne") {
        if (!requireNamespace("Rtsne", quietly = TRUE)) {
          warning("Package 'Rtsne' required for t-SNE. Using PCA instead.")
          embeddings <- corrected$corrected[, 1:2]
        } else {
          tsne_res <- Rtsne::Rtsne(corrected$corrected, dims = 2, perplexity = 30)$Y
          embeddings <- tsne_res
          SingleCellExperiment::reducedDim(integrated_data, "TSNE") <- tsne_res
        }
      } else {
        embeddings <- corrected$corrected[, 1:2]
      }

      # Store results
      results$integration_results <- integrated_data
      results$integration_embeddings <- embeddings

      # Create plots if ggplot2 is available
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        # Combine metadata from all modalities
        metadata <- do.call(rbind, lapply(1:n_modalities, function(i) {
          md <- as.data.frame(SingleCellExperiment::colData(sce_list[[i]]))
          md$modality <- modality_names[i]
          return(md)
        }))

        # Create dimension reduction plot
        if (cell_type_var %in% colnames(metadata)) {
          plot_data <- data.frame(
            x = embeddings[, 1],
            y = embeddings[, 2],
            cell_type = metadata[[cell_type_var]]
          )

          dim_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = cell_type)) +
            ggplot2::geom_point(size = 1, alpha = 0.7) +
            ggplot2::theme_minimal() +
            ggplot2::labs(
              title = paste0("Integrated ", dim_reduction, " Plot"),
              x = "Dimension 1",
              y = "Dimension 2",
              color = "Cell Type"
            )

          results$plots$dim_plot <- dim_plot
        }
      }
    } else {
      stop("All objects must be of the same type (either all Seurat or all SingleCellExperiment)")
    }

  } else if (integration_method == "mnn") {
    # Mutual Nearest Neighbors (MNN) integration

    if (!requireNamespace("batchelor", quietly = TRUE)) {
      stop("Package 'batchelor' required for MNN integration")
    }

    if (all(object_types == "Seurat")) {
      # Convert Seurat objects to matrices for batchelor
      expression_matrices <- lapply(data_list, function(seurat_obj) {
        return(as.matrix(Seurat::GetAssayData(seurat_obj, slot = "data")))
      })

      # Perform MNN batch correction
      message("Performing MNN integration...")
      mnn_results <- batchelor::fastMNN(
        expression_matrices,
        d = n_dims,
        k = 20,
        auto.merge = TRUE
      )

      # Create a Seurat object from corrected data
      integrated_data <- Seurat::CreateSeuratObject(
        counts = expression_matrices[[1]],
        meta.data = data_list[[1]]@meta.data
      )

      # Add corrected dimensions
      integrated_data[["mnn"]] <- Seurat::CreateDimReducObject(
        embeddings = mnn_results$corrected,
        key = "MNN_"
      )

      # Run dimension reduction on corrected data
      if (dim_reduction == "umap") {
        integrated_data <- Seurat::RunUMAP(integrated_data, dims = 1:n_dims, reduction = "mnn")
        embeddings <- Seurat::Embeddings(integrated_data, reduction = "umap")
      } else if (dim_reduction == "tsne") {
        integrated_data <- Seurat::RunTSNE(integrated_data, dims = 1:n_dims, reduction = "mnn")
        embeddings <- Seurat::Embeddings(integrated_data, reduction = "tsne")
      } else {
        embeddings <- mnn_results$corrected[, 1:2]
      }

    } else if (all(object_types == "SCE")) {
      # Perform MNN batch correction directly on SCE objects
      message("Performing MNN integration with SingleCellExperiment objects...")

      sce_list <- data_list

      # Perform MNN batch correction
      corrected <- batchelor::fastMNN(
        sce_list,
        d = n_dims,
        k = 20,
        auto.merge = TRUE
      )

      # Store corrected dimensions in a new SingleCellExperiment
      integrated_data <- corrected

      # Run dimension reduction on corrected data
      if (dim_reduction == "umap") {
        if (!requireNamespace("scater", quietly = TRUE)) {
          warning("Package 'scater' required for UMAP. Using PCA instead.")
          embeddings <- SingleCellExperiment::reducedDim(integrated_data, "corrected")[, 1:2]
        } else {
          integrated_data <- scater::runUMAP(integrated_data, dimred = "corrected")
          embeddings <- SingleCellExperiment::reducedDim(integrated_data, "UMAP")
        }
      } else if (dim_reduction == "tsne") {
        if (!requireNamespace("scater", quietly = TRUE)) {
          warning("Package 'scater' required for t-SNE. Using PCA instead.")
          embeddings <- SingleCellExperiment::reducedDim(integrated_data, "corrected")[, 1:2]
        } else {
          integrated_data <- scater::runTSNE(integrated_data, dimred = "corrected")
          embeddings <- SingleCellExperiment::reducedDim(integrated_data, "TSNE")
        }
      } else {
        embeddings <- SingleCellExperiment::reducedDim(integrated_data, "corrected")[, 1:2]
      }
    } else {
      stop("All objects must be of the same type (either all Seurat or all SingleCellExperiment)")
    }

    # Store results
    results$integration_results <- integrated_data
    results$integration_embeddings <- embeddings

  } else if (integration_method == "harmony") {
    # Harmony integration

    if (!requireNamespace("harmony", quietly = TRUE)) {
      stop("Package 'harmony' required for Harmony integration")
    }

    if (all(object_types == "Seurat")) {
      # Merge Seurat objects
      merged_object <- data_list[[1]]

      # Add modality information to metadata
      merged_object$modality <- modality_names[1]

      for (i in 2:n_modalities) {
        # Add modality information to each object
        data_list[[i]]$modality <- modality_names[i]

        # Merge with the first object
        merged_object <- merge(merged_object, data_list[[i]])
      }

      # Process merged object
      merged_object <- Seurat::FindVariableFeatures(merged_object,
                                                    selection.method = "vst",
                                                    nfeatures = n_features)
      merged_object <- Seurat::ScaleData(merged_object)
      merged_object <- Seurat::RunPCA(merged_object, npcs = n_dims)

      # Run Harmony integration
      message("Performing Harmony integration...")
      harmony_embeddings <- harmony::RunHarmony(
        merged_object,
        group.by.vars = "modality",
        reduction = "pca",
        dims.use = 1:n_dims,
        max.iter.harmony = 20
      )

      # Add Harmony embeddings to the object
      merged_object[["harmony"]] <- harmony_embeddings

      # Run dimension reduction on Harmony embeddings
      if (dim_reduction == "umap") {
        merged_object <- Seurat::RunUMAP(merged_object, reduction = "harmony", dims = 1:n_dims)
        embeddings <- Seurat::Embeddings(merged_object, reduction = "umap")
      } else if (dim_reduction == "tsne") {
        merged_object <- Seurat::RunTSNE(merged_object, reduction = "harmony", dims = 1:n_dims)
        embeddings <- Seurat::Embeddings(merged_object, reduction = "tsne")
      } else {
        embeddings <- Seurat::Embeddings(merged_object, reduction = "harmony")[, 1:2]
      }

      integrated_data <- merged_object

    } else if (all(object_types == "SCE")) {
      # Combine SCE objects
      combined_sce <- SingleCellExperiment::cbind(data_list[[1]], data_list[[2]])

      # Add modality information to colData
      combined_sce$modality <- rep(modality_names[1:2], c(ncol(data_list[[1]]), ncol(data_list[[2]])))

      if (n_modalities > 2) {
        for (i in 3:n_modalities) {
          # Add additional SCE objects
          temp_sce <- combined_sce
          combined_sce <- SingleCellExperiment::cbind(temp_sce, data_list[[i]])
          combined_sce$modality <- c(temp_sce$modality, rep(modality_names[i], ncol(data_list[[i]])))
        }
      }

      # Run PCA on combined data
      if (!requireNamespace("scater", quietly = TRUE)) {
        stop("Package 'scater' required for processing SingleCellExperiment objects")
      }

      combined_sce <- scater::runPCA(combined_sce, ncomponents = n_dims)

      # Run Harmony on PCA embeddings
      message("Performing Harmony integration with SingleCellExperiment objects...")
      pca_mat <- SingleCellExperiment::reducedDim(combined_sce, "PCA")

      harmony_out <- harmony::HarmonyMatrix(
        pca_mat,
        combined_sce$modality,
        do_pca = FALSE,
        npcs = n_dims
      )

      # Store Harmony results
      SingleCellExperiment::reducedDim(combined_sce, "harmony") <- harmony_out

      # Run dimension reduction on Harmony embeddings
      if (dim_reduction == "umap") {
        combined_sce <- scater::runUMAP(combined_sce, dimred = "harmony")
        embeddings <- SingleCellExperiment::reducedDim(combined_sce, "UMAP")
      } else if (dim_reduction == "tsne") {
        combined_sce <- scater::runTSNE(combined_sce, dimred = "harmony")
        embeddings <- SingleCellExperiment::reducedDim(combined_sce, "TSNE")
      } else {
        embeddings <- harmony_out[, 1:2]
      }

      integrated_data <- combined_sce

    } else {
      stop("All objects must be of the same type (either all Seurat or all SingleCellExperiment)")
    }

    # Store results
    results$integration_results <- integrated_data
    results$integration_embeddings <- embeddings

  } else if (integration_method == "seurat") {
    # Seurat integration (only available for Seurat objects)

    if (!all(object_types == "Seurat")) {
      stop("Seurat integration method is only available for Seurat objects")
    }

    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' required for Seurat integration")
    }

    # Perform standard Seurat integration workflow
    message("Performing Seurat integration...")

    # Normalize and identify variable features for each dataset
    data_list <- lapply(data_list, function(x) {
      x <- Seurat::NormalizeData(x)
      x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = n_features)
      return(x)
    })

    # Select features for integration
    features <- Seurat::SelectIntegrationFeatures(data_list, nfeatures = n_features)

    # Find integration anchors
    anchors <- Seurat::FindIntegrationAnchors(
      object.list = data_list,
      anchor.features = features,
      normalization.method = "LogNormalize"
    )

    # Integrate data
    integrated_data <- Seurat::IntegrateData(anchors)

    # Switch to integrated assay and scale data
    Seurat::DefaultAssay(integrated_data) <- "integrated"
    integrated_data <- Seurat::ScaleData(integrated_data)

    # Run dimension reduction
    integrated_data <- Seurat::RunPCA(integrated_data, npcs = n_dims)

    if (dim_reduction == "umap") {
      integrated_data <- Seurat::RunUMAP(integrated_data, dims = 1:n_dims)
      embeddings <- Seurat::Embeddings(integrated_data, reduction = "umap")
    } else if (dim_reduction == "tsne") {
      integrated_data <- Seurat::RunTSNE(integrated_data, dims = 1:n_dims)
      embeddings <- Seurat::Embeddings(integrated_data, reduction = "tsne")
    } else {
      embeddings <- Seurat::Embeddings(integrated_data, reduction = "pca")[, 1:2]
    }

    # Store results
    results$integration_results <- integrated_data
    results$integration_embeddings <- embeddings

  } else {
    stop(paste0("Unsupported integration method: ", integration_method))
  }

  # Generate visualization if embeddings are available
  if (!is.null(results$integration_embeddings) && requireNamespace("ggplot2", quietly = TRUE)) {
    # Extract cell metadata
    if (inherits(results$integration_results, "Seurat")) {
      metadata <- results$integration_results@meta.data
    } else if (inherits(results$integration_results, "SingleCellExperiment")) {
      metadata <- as.data.frame(SingleCellExperiment::colData(results$integration_results))
    }

    # Create embedding visualization
    embeddings_df <- data.frame(
      Dim1 = results$integration_embeddings[, 1],
      Dim2 = results$integration_embeddings[, 2]
    )

    # Add modality and cell type information if available
    if ("modality" %in% colnames(metadata)) {
      embeddings_df$Modality <- metadata$modality
    }

    if (cell_type_var %in% colnames(metadata)) {
      embeddings_df$CellType <- metadata[[cell_type_var]]
    }

    # Create plots
    if ("Modality" %in% colnames(embeddings_df)) {
      modality_plot <- ggplot2::ggplot(embeddings_df, ggplot2::aes(x = Dim1, y = Dim2, color = Modality)) +
        ggplot2::geom_point(size = 1, alpha = 0.7) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste0("Integration by ", integration_method, " - Colored by Modality"),
          x = "Dimension 1",
          y = "Dimension 2"
        )

      results$plots$modality_plot <- modality_plot
    }

    if ("CellType" %in% colnames(embeddings_df)) {
      celltype_plot <- ggplot2::ggplot(embeddings_df, ggplot2::aes(x = Dim1, y = Dim2, color = CellType)) +
        ggplot2::geom_point(size = 1, alpha = 0.7) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste0("Integration by ", integration_method, " - Colored by Cell Type"),
          x = "Dimension 1",
          y = "Dimension 2"
        )

      results$plots$celltype_plot <- celltype_plot
    }
  }

  return(results)
}

#' Calculate Agreement Between Modalities
#'
#' Calculates the agreement or concordance between different modalities after integration.
#'
#' @param integration_results Output from integrateMultiModalData function
#' @param features_to_compare Features to compare across modalities (default: variable features)
#' @param correlation_method Method for correlation calculation: "pearson", "spearman", or "kendall" (default: "spearman")
#' @param n_neighbors Number of neighbors to consider for KNN-based metrics (default: 15)
#' @param cell_subset Optional vector of cell indices or names to subset for analysis
#'
#' @return A list containing agreement metrics and visualizations
#' @export
calculateModalityAgreement <- function(integration_results,
                                       features_to_compare = NULL,
                                       correlation_method = "spearman",
                                       n_neighbors = 15,
                                       cell_subset = NULL) {

  # Extract integrated data and embeddings
  integrated_data <- integration_results$integration_results
  embeddings <- integration_results$integration_embeddings

  # Check if input is valid
  if (is.null(integrated_data)) {
    stop("No integration results found")
  }

  # Determine object type
  is_seurat <- inherits(integrated_data, "Seurat")
  is_sce <- inherits(integrated_data, "SingleCellExperiment")

  if (!is_seurat && !is_sce) {
    stop("Integrated data must be a Seurat or SingleCellExperiment object")
  }

  # Apply cell subset if provided
  if (!is.null(cell_subset)) {
    if (is_seurat) {
      integrated_data <- subset(integrated_data, cells = cell_subset)
    } else {
      integrated_data <- integrated_data[, cell_subset]
    }
    embeddings <- embeddings[cell_subset, ]
  }

  # Initialize results
  agreement_results <- list(
    metadata = list(
      correlation_method = correlation_method,
      n_neighbors = n_neighbors
    ),
    feature_correlations = NULL,
    neighborhood_stability = NULL,
    modality_dist_correlations = NULL,
    plots = list()
  )

  # Extract modality information
  if (is_seurat) {
    metadata <- integrated_data@meta.data
    modalities <- unique(metadata$modality)
  } else {
    metadata <- as.data.frame(SingleCellExperiment::colData(integrated_data))
    modalities <- unique(metadata$modality)
  }

  # Compute feature-level correlations
  if (!is.null(features_to_compare)) {
    # Extract expression matrices for each modality
    modality_expressions <- list()

    for (mod in modalities) {
      if (is_seurat) {
        cells <- rownames(metadata)[metadata$modality == mod]
        if (length(cells) > 0) {
          # Get expression data
          if ("integrated" %in% Seurat::Assays(integrated_data)) {
            expr <- as.matrix(Seurat::GetAssayData(integrated_data, slot = "data", assay = "integrated"))
          } else {
            expr <- as.matrix(Seurat::GetAssayData(integrated_data, slot = "data"))
          }
          modality_expressions[[mod]] <- expr[features_to_compare, cells]
        }
      } else {
        cells <- which(metadata$modality == mod)
        if (length(cells) > 0) {
          # Get expression data
          if ("logcounts" %in% SingleCellExperiment::assayNames(integrated_data)) {
            expr <- SingleCellExperiment::logcounts(integrated_data)
          } else {
            expr <- SingleCellExperiment::counts(integrated_data)
          }
          modality_expressions[[mod]] <- expr[features_to_compare, cells]
        }
      }
    }

    # Calculate correlation between modalities for each feature
    feature_cors <- matrix(NA, nrow = length(features_to_compare),
                           ncol = length(modalities) * (length(modalities) - 1) / 2)
    rownames(feature_cors) <- features_to_compare

    # Create column names for modality pairs
    mod_pairs <- combn(modalities, 2)
    colnames(feature_cors) <- apply(mod_pairs, 2, paste, collapse = "_vs_")

    # Compute correlations
    for (i in 1:ncol(mod_pairs)) {
      mod1 <- mod_pairs[1, i]
      mod2 <- mod_pairs[2, i]

      # Skip if data not available
      if (!mod1 %in% names(modality_expressions) || !mod2 %in% names(modality_expressions)) {
        next
      }

      expr1 <- modality_expressions[[mod1]]
      expr2 <- modality_expressions[[mod2]]

      # Compute correlation for each feature
      for (j in 1:length(features_to_compare)) {
        feature <- features_to_compare[j]
        if (feature %in% rownames(expr1) && feature %in% rownames(expr2)) {
          feature_cors[j, i] <- cor(t(expr1[feature, , drop = FALSE]),
                                    t(expr2[feature, , drop = FALSE]),
                                    method = correlation_method)
        }
      }
    }

    # Store feature correlations
    agreement_results$feature_correlations <- feature_cors

    # Create heatmap visualization
    if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("reshape2", quietly = TRUE)) {
      # Reshape data for heatmap
      feature_cors_long <- reshape2::melt(feature_cors, varnames = c("Feature", "Comparison"),
                                          value.name = "Correlation")

      # Create heatmap
      feature_heatmap <- ggplot2::ggplot(feature_cors_long, ggplot2::aes(x = Comparison, y = Feature,
                                                                         fill = Correlation)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                                      limits = c(-1, 1), na.value = "grey90") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(title = "Feature Correlation Between Modalities",
                      x = "Modality Comparison", y = "Feature")

      agreement_results$plots$feature_heatmap <- feature_heatmap
    }
  }

  # Compute neighborhood stability
  # This measures whether cells have similar neighbors in different modalities
  if (!is.null(embeddings) && requireNamespace("FNN", quietly = TRUE)) {
    # Get KNN for each cell in the integrated space
    knn_integrated <- FNN::get.knn(embeddings, k = n_neighbors)

    # Get cells from each modality
    mod_cells <- list()
    for (mod in modalities) {
      if (is_seurat) {
        mod_cells[[mod]] <- rownames(metadata)[metadata$modality == mod]
      } else {
        mod_cells[[mod]] <- rownames(metadata)[metadata$modality == mod]
      }
    }

    # Initialize neighborhood stability matrix
    neighborhood_stability <- matrix(NA, nrow = nrow(embeddings), ncol = length(modalities))
    rownames(neighborhood_stability) <- rownames(embeddings)
    colnames(neighborhood_stability) <- modalities

    # For each cell, compute stability metric
    for (i in 1:nrow(embeddings)) {
      # Get KNN indices for this cell
      nn_indices <- knn_integrated$nn.index[i, ]

      # For each modality, compute overlap with KNN
      for (m in 1:length(modalities)) {
        mod <- modalities[m]
        # Skip if no cells for this modality
        if (length(mod_cells[[mod]]) == 0) {
          next
        }

        # Count neighbors from this modality
        neighbors_from_mod <- sum(rownames(embeddings)[nn_indices] %in% mod_cells[[mod]])

        # Compute stability as proportion of neighbors from this modality
        neighborhood_stability[i, m] <- neighbors_from_mod / n_neighbors
      }
    }

    # Store neighborhood stability
    agreement_results$neighborhood_stability <- neighborhood_stability

    # Create neighborhood stability visualization
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      # Compute average stability for each modality
      avg_stability <- colMeans(neighborhood_stability, na.rm = TRUE)

      # Create barplot
      stability_plot <- ggplot2::ggplot(data.frame(Modality = modalities,
                                                   Stability = avg_stability),
                                        ggplot2::aes(x = Modality, y = Stability, fill = Modality)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Neighborhood Stability by Modality",
                      x = "Modality", y = "Average Stability")

      agreement_results$plots$stability_plot <- stability_plot
    }
  }

  # Calculate distance correlation between modalities
  # This compares distance matrices in each modality's native space
  if (length(modalities) > 1) {
    # Initialize distance correlation matrix
    dist_correlation <- matrix(NA, nrow = length(modalities), ncol = length(modalities))
    rownames(dist_correlation) <- modalities
    colnames(dist_correlation) <- modalities

    # For each pair of modalities, compute distance correlation
    for (i in 1:length(modalities)) {
      for (j in 1:length(modalities)) {
        mod_i <- modalities[i]
        mod_j <- modalities[j]

        # Skip diagonal or if data not available
        if (i == j || !mod_i %in% names(modality_expressions) || !mod_j %in% names(modality_expressions)) {
          if (i == j) {
            dist_correlation[i, j] <- 1  # Set diagonal to 1
          }
          next
        }

        # Compute distance matrices
        if (requireNamespace("stats", quietly = TRUE)) {
          # Get common features
          common_features <- intersect(rownames(modality_expressions[[mod_i]]),
                                       rownames(modality_expressions[[mod_j]]))

          if (length(common_features) == 0) {
            next
          }

          # Subset to common features
          expr_i <- modality_expressions[[mod_i]][common_features, ]
          expr_j <- modality_expressions[[mod_j]][common_features, ]

          # Compute distance matrices (cell-cell distances)
          dist_i <- stats::dist(t(expr_i))
          dist_j <- stats::dist(t(expr_j))

          # Convert to matrices for correlation
          dist_mat_i <- as.matrix(dist_i)
          dist_mat_j <- as.matrix(dist_j)

          # Compute correlation between distance matrices
          dist_correlation[i, j] <- cor(dist_mat_i[lower.tri(dist_mat_i)],
                                        dist_mat_j[lower.tri(dist_mat_j)],
                                        method = correlation_method)
        }
      }
    }

    # Store distance correlation
    agreement_results$modality_dist_correlations <- dist_correlation

    # Create distance correlation heatmap
    if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("reshape2", quietly = TRUE)) {
      # Reshape data for heatmap
      dist_cor_long <- reshape2::melt(dist_correlation, varnames = c("Modality1", "Modality2"),
                                      value.name = "Correlation")

      # Create heatmap
      dist_heatmap <- ggplot2::ggplot(dist_cor_long, ggplot2::aes(x = Modality1, y = Modality2,
                                                                  fill = Correlation)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5,
                                      limits = c(0, 1), na.value = "grey90") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Distance Matrix Correlation Between Modalities",
                      x = "Modality 1", y = "Modality 2")

      agreement_results$plots$dist_heatmap <- dist_heatmap
    }
  }

  return(agreement_results)
}

#' Find Conserved Features Across Modalities
#'
#' Identifies features that show consistent patterns across different data modalities.
#'
#' @param integration_results Output from integrateMultiModalData function
#' @param cluster_var Name of the column containing cluster assignments
#' @param min_cells_per_cluster Minimum number of cells required per cluster (default: 10)
#' @param p_val_threshold P-value threshold for significance (default: 0.05)
#' @param min_pct Minimum percentage of cells expressing the feature (default: 0.1)
#' @param min_diff_pct Minimum percentage difference between the two groups of cells (default: 0.1)
#' @param only_positive If TRUE, only return positively conserved features (default: TRUE)
#'
#' @return A list of conserved features for each cluster across modalities
#' @export
findConservedFeatures <- function(integration_results,
                                  cluster_var,
                                  min_cells_per_cluster = 10,
                                  p_val_threshold = 0.05,
                                  min_pct = 0.1,
                                  min_diff_pct = 0.1,
                                  only_positive = TRUE) {

  # Extract integrated data
  integrated_data <- integration_results$integration_results

  # Check if input is valid
  if (is.null(integrated_data)) {
    stop("No integration results found")
  }

  # Determine object type
  is_seurat <- inherits(integrated_data, "Seurat")
  is_sce <- inherits(integrated_data, "SingleCellExperiment")

  if (!is_seurat && !is_sce) {
    stop("Integrated data must be a Seurat or SingleCellExperiment object")
  }

  # Extract metadata
  if (is_seurat) {
    metadata <- integrated_data@meta.data
  } else {
    metadata <- as.data.frame(SingleCellExperiment::colData(integrated_data))
  }

  # Check if required columns exist in metadata
  if (!cluster_var %in% colnames(metadata)) {
    stop(paste0("Column '", cluster_var, "' not found in metadata"))
  }

  if (!"modality" %in% colnames(metadata)) {
    stop("Column 'modality' not found in metadata")
  }

  # Extract modalities and clusters
  modalities <- unique(metadata$modality)
  clusters <- unique(metadata[[cluster_var]])

  if (length(modalities) < 2) {
    stop("At least two modalities are required for conserved feature analysis")
  }

  # Initialize results
  conserved_features <- list()

  # Process each cluster
  for (cluster in clusters) {
    # Get cells in this cluster
    cluster_cells <- rownames(metadata)[metadata[[cluster_var]] == cluster]

    # Skip if too few cells
    if (length(cluster_cells) < min_cells_per_cluster) {
      message(paste0("Skipping cluster ", cluster, ": insufficient cells (",
                     length(cluster_cells), " < ", min_cells_per_cluster, ")"))
      next
    }

    # Run analysis based on object type
    if (is_seurat) {
      # Use Seurat's FindConservedMarkers
      if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' required for finding conserved features")
      }

      # Set identity to cluster variable
      Seurat::Idents(integrated_data) <- integrated_data@meta.data[[cluster_var]]

      # Find conserved markers
      conserved_markers <- Seurat::FindConservedMarkers(
        integrated_data,
        ident.1 = cluster,
        grouping.var = "modality",
        only.pos = only_positive,
        min.pct = min_pct,
        min.diff.pct = min_diff_pct,
        logfc.threshold = 0.25
      )

      # Filter by p-value threshold
      if (!is.null(conserved_markers) && nrow(conserved_markers) > 0) {
        # Check if any feature passes the threshold in all modalities
        max_p_val_cols <- grep("_p_val$", colnames(conserved_markers))
        if (length(max_p_val_cols) > 0) {
          valid_features <- apply(conserved_markers[, max_p_val_cols, drop = FALSE], 1,
                                  function(x) all(x < p_val_threshold))
          conserved_markers <- conserved_markers[valid_features, ]
        }
      }

      # Store results
      if (!is.null(conserved_markers) && nrow(conserved_markers) > 0) {
        conserved_markers$gene <- rownames(conserved_markers)
        conserved_features[[cluster]] <- conserved_markers
      }

    } else {
      # For SingleCellExperiment, implement custom conserved feature detection
      # First, get expression data
      if ("logcounts" %in% SingleCellExperiment::assayNames(integrated_data)) {
        expr_data <- SingleCellExperiment::logcounts(integrated_data)
      } else {
        expr_data <- SingleCellExperiment::counts(integrated_data)
      }

      # Initialize results for this cluster
      cluster_results <- data.frame(gene = rownames(expr_data))

      # Process each modality
      for (mod in modalities) {
        # Get cells from this modality
        mod_cells <- rownames(metadata)[metadata$modality == mod]

        # Skip if too few cells
        if (length(mod_cells) < 5) {
          next
        }

        # Split cells into two groups: in cluster vs not in cluster
        in_cluster_cells <- intersect(cluster_cells, mod_cells)
        out_cluster_cells <- setdiff(mod_cells, in_cluster_cells)

        # Skip if too few cells in either group
        if (length(in_cluster_cells) < 3 || length(out_cluster_cells) < 3) {
          next
        }

        # Calculate percent expressed
        pct_in <- rowMeans(expr_data[, in_cluster_cells, drop = FALSE] > 0)
        pct_out <- rowMeans(expr_data[, out_cluster_cells, drop = FALSE] > 0)

        # Calculate average expression
        avg_in <- rowMeans(expr_data[, in_cluster_cells, drop = FALSE])
        avg_out <- rowMeans(expr_data[, out_cluster_cells, drop = FALSE])

        # Calculate log fold change
        logfc <- avg_in - avg_out

        # Perform statistical test (Wilcoxon)
        p_vals <- numeric(length(rownames(expr_data)))
        names(p_vals) <- rownames(expr_data)

        for (i in 1:length(rownames(expr_data))) {
          gene <- rownames(expr_data)[i]
          expr_in <- expr_data[gene, in_cluster_cells]
          expr_out <- expr_data[gene, out_cluster_cells]

          # Skip if too little expression
          if (mean(expr_in > 0) < min_pct || abs(mean(expr_in > 0) - mean(expr_out > 0)) < min_diff_pct) {
            p_vals[i] <- 1
            next
          }

          # Wilcoxon test
          tryCatch({
            test_result <- stats::wilcox.test(expr_in, expr_out)
            p_vals[i] <- test_result$p.value
          }, error = function(e) {
            p_vals[i] <- 1
          })
        }

        # Adjust p-values
        p_adj <- stats::p.adjust(p_vals, method = "BH")

        # Store results for this modality
        cluster_results[[paste0(mod, "_logfc")]] <- logfc
        cluster_results[[paste0(mod, "_pct.1")]] <- pct_in
        cluster_results[[paste0(mod, "_pct.2")]] <- pct_out
        cluster_results[[paste0(mod, "_p_val")]] <- p_vals
        cluster_results[[paste0(mod, "_p_val_adj")]] <- p_adj
      }

      # Filter for conserved features across modalities
      p_val_cols <- grep("_p_val_adj$", colnames(cluster_results))
      logfc_cols <- grep("_logfc$", colnames(cluster_results))

      if (length(p_val_cols) > 0 && length(logfc_cols) > 0) {
        # Check all modalities pass p-value threshold
        valid_p_val <- apply(cluster_results[, p_val_cols, drop = FALSE], 1,
                             function(x) all(x < p_val_threshold))

        # Check all modalities have positive log fold change (if only_positive)
        if (only_positive) {
          valid_logfc <- apply(cluster_results[, logfc_cols, drop = FALSE], 1,
                               function(x) all(x > 0))
        } else {
          valid_logfc <- apply(cluster_results[, logfc_cols, drop = FALSE], 1,
                               function(x) all(sign(x) == sign(x[1])))
        }

        # Combine criteria
        valid_features <- valid_p_val & valid_logfc

        # Filter results
        if (sum(valid_features) > 0) {
          conserved_features[[cluster]] <- cluster_results[valid_features, ]
        }
      }
    }
  }

  # Return results
  return(conserved_features)
}

#' Project New Data onto Integrated Space
#'
#' Projects data from a new experiment onto the existing integrated space.
#'
#' @param integration_results Output from integrateMultiModalData function
#' @param new_data New data object (must be same class as original data)
#' @param reference_modality Name of the modality to use as reference
#' @param variable_features Optional vector of features to use for projection
#' @param scale_new_data Logical, whether to scale new data (default: TRUE)
#'
#' @return A list containing projected data and visualization
#' @export
projectNewData <- function(integration_results,
                           new_data,
                           reference_modality,
                           variable_features = NULL,
                           scale_new_data = TRUE) {

  # Extract integrated data and embeddings
  integrated_data <- integration_results$integration_results
  embeddings <- integration_results$integration_embeddings

  # Check if input is valid
  if (is.null(integrated_data)) {
    stop("No integration results found")
  }

  # Determine object type and validate new data
  is_seurat <- inherits(integrated_data, "Seurat")
  is_sce <- inherits(integrated_data, "SingleCellExperiment")

  if (!is_seurat && !is_sce) {
    stop("Integrated data must be a Seurat or SingleCellExperiment object")
  }

  if (is_seurat && !inherits(new_data, "Seurat")) {
    stop("New data must be the same class as original data (Seurat)")
  }

  if (is_sce && !inherits(new_data, "SingleCellExperiment")) {
    stop("New data must be the same class as original data (SingleCellExperiment)")
  }

  # Get modality information
  if (is_seurat) {
    metadata <- integrated_data@meta.data
    modalities <- unique(metadata$modality)
  } else {
    metadata <- as.data.frame(SingleCellExperiment::colData(integrated_data))
    modalities <- unique(metadata$modality)
  }

  # Check if reference modality exists
  if (!reference_modality %in% modalities) {
    stop(paste0("Reference modality '", reference_modality, "' not found in integration results"))
  }

  # Initialize results
  projection_results <- list(
    metadata = list(
      reference_modality = reference_modality,
      variable_features = variable_features
    ),
    projected_data = NULL,
    projected_embeddings = NULL,
    plots = list()
  )

  # Handle variable features
  if (is.null(variable_features)) {
    if (is_seurat) {
      variable_features <- Seurat::VariableFeatures(integrated_data)
    } else {
      # For SCE, try to get variable features from metadata or compute them
      if ("highly_variable" %in% names(SingleCellExperiment::rowData(integrated_data))) {
        variable_features <- rownames(integrated_data)[SingleCellExperiment::rowData(integrated_data)$highly_variable]
      } else if (requireNamespace("scran", quietly = TRUE)) {
        message("Computing variable features...")
        var_features <- scran::modelGeneVar(integrated_data)
        variable_features <- rownames(var_features)[order(var_features$bio, decreasing = TRUE)[1:2000]]
      } else {
        # Default to top 2000 genes by variance
        expr <- SingleCellExperiment::logcounts(integrated_data)
        gene_vars <- apply(expr, 1, stats::var)
        variable_features <- names(sort(gene_vars, decreasing = TRUE)[1:2000])
      }
    }
  }

  # Project new data based on object type
  if (is_seurat) {
    # Use Seurat's projection functions
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' required for data projection")
    }

    # Process based on integration method
    integration_method <- integration_results$metadata$integration_method

    if (integration_method == "harmony") {
      # For Harmony integration
      if (!requireNamespace("harmony", quietly = TRUE)) {
        stop("Package 'harmony' required for projection with Harmony integration")
      }

      # Process new data
      new_data <- Seurat::FindVariableFeatures(new_data, selection.method = "vst",
                                               nfeatures = length(variable_features))
      new_data <- Seurat::ScaleData(new_data, features = variable_features)
      new_data <- Seurat::RunPCA(new_data, features = variable_features, npcs = 50)

      # Get harmony embeddings
      harmony_embedding <- Seurat::Embeddings(integrated_data, reduction = "harmony")

      # Get reference cells
      ref_cells <- rownames(metadata)[metadata$modality == reference_modality]

      # Project new data to harmony space
      projected_harmony <- projectToLowDim(
        new_data = Seurat::Embeddings(new_data, reduction = "pca"),
        ref_data = Seurat::Embeddings(integrated_data, reduction = "pca")[ref_cells, ],
        ref_embedding = harmony_embedding[ref_cells, ]
      )

      # Create new reduction
      new_data[["harmony"]] <- Seurat::CreateDimReducObject(
        embeddings = projected_harmony,
        key = "harmony_"
      )

      # Run UMAP on projected data
      if ("umap" %in% names(integrated_data@reductions)) {
        new_data <- Seurat::RunUMAP(new_data, reduction = "harmony",
                                    dims = 1:ncol(projected_harmony),
                                    reduction.name = "umap")
      }

      # Store projected data
      projection_results$projected_data <- new_data
      projection_results$projected_embeddings <- Seurat::Embeddings(new_data, reduction = "harmony")

    } else if (integration_method %in% c("canonical_correlation", "seurat")) {
      # For Seurat integration methods
      # Find anchors between reference and new data
      reference_data <- subset(integrated_data, cells = rownames(metadata)[metadata$modality == reference_modality])

      # Find anchors
      transfer_anchors <- Seurat::FindTransferAnchors(
        reference = reference_data,
        query = new_data,
        dims = 1:30,
        reference.reduction = "pca"
      )

      # Project data
      projected_data <- Seurat::MapQuery(
        anchorset = transfer_anchors,
        reference = reference_data,
        query = new_data,
        refdata = list(cell_type = reference_data@meta.data[[colnames(reference_data@meta.data)[1]]]),
        reference.reduction = "pca",
        reduction.model = "umap"
      )

      # Store projected data
      projection_results$projected_data <- projected_data
      projection_results$projected_embeddings <- Seurat::Embeddings(projected_data, reduction = "ref.pca")

    } else if (integration_method == "mnn") {
      # For MNN integration
      if (!requireNamespace("batchelor", quietly = TRUE)) {
        stop("Package 'batchelor' required for projection with MNN integration")
      }

      # Process new data
      new_data <- Seurat::FindVariableFeatures(new_data, selection.method = "vst",
                                               nfeatures = length(variable_features))
      new_data <- Seurat::ScaleData(new_data, features = variable_features)

      # Get reference cells
      ref_cells <- rownames(metadata)[metadata$modality == reference_modality]
      ref_data <- subset(integrated_data, cells = ref_cells)

      # Extract expression matrices
      new_expr <- as.matrix(Seurat::GetAssayData(new_data, slot = "data"))
      ref_expr <- as.matrix(Seurat::GetAssayData(ref_data, slot = "data"))

      # Find common genes
      common_genes <- intersect(rownames(new_expr), rownames(ref_expr))

      # Project using batchelor
      projected_mnn <- batchelor::fastMNN(
        new_expr[common_genes, ],
        ref_expr[common_genes, ],
        d = 30
      )

      # Create new reduction
      new_data[["mnn"]] <- Seurat::CreateDimReducObject(
        embeddings = projected_mnn$corrected,
        key = "mnn_"
      )

      # Run UMAP on projected data
      if ("umap" %in% names(integrated_data@reductions)) {
        new_data <- Seurat::RunUMAP(new_data, reduction = "mnn",
                                    dims = 1:ncol(projected_mnn$corrected),
                                    reduction.name = "umap")
      }

      # Store projected data
      projection_results$projected_data <- new_data
      projection_results$projected_embeddings <- projected_mnn$corrected
    }

  } else {
    # For SingleCellExperiment objects
    # Get reference cells
    ref_cells <- rownames(metadata)[metadata$modality == reference_modality]
    ref_data <- integrated_data[, ref_cells]

    # Process new data
    if (scale_new_data) {
      # Scale new data
      if (requireNamespace("scater", quietly = TRUE)) {
        new_data <- scater::logNormCounts(new_data)
      }
    }

    # Get expression data
    if ("logcounts" %in% SingleCellExperiment::assayNames(new_data)) {
      new_expr <- SingleCellExperiment::logcounts(new_data)
    } else {
      new_expr <- SingleCellExperiment::counts(new_data)
    }

    if ("logcounts" %in% SingleCellExperiment::assayNames(ref_data)) {
      ref_expr <- SingleCellExperiment::logcounts(ref_data)
    } else {
      ref_expr <- SingleCellExperiment::counts(ref_data)
    }

    # Find common genes
    common_genes <- intersect(rownames(new_expr), rownames(ref_expr))

    # Project data based on integration method
    integration_method <- integration_results$metadata$integration_method

    if (integration_method == "harmony") {
      # For Harmony integration
      if (!requireNamespace("harmony", quietly = TRUE)) {
        stop("Package 'harmony' required for projection with Harmony integration")
      }

      # Run PCA on new data
      if (requireNamespace("scater", quietly = TRUE)) {
        new_data <- scater::runPCA(new_data, ncomponents = 50)
      }

      # Get PCA and harmony embeddings
      pca_embedding <- SingleCellExperiment::reducedDim(ref_data, "PCA")
      harmony_embedding <- SingleCellExperiment::reducedDim(ref_data, "harmony")

      # Project new data to harmony space
      projected_harmony <- projectToLowDim(
        new_data = SingleCellExperiment::reducedDim(new_data, "PCA"),
        ref_data = pca_embedding,
        ref_embedding = harmony_embedding
      )

      # Add projected dimensions
      SingleCellExperiment::reducedDim(new_data, "harmony") <- projected_harmony

      # Run UMAP on projected data
      if ("UMAP" %in% SingleCellExperiment::reducedDimNames(ref_data)) {
        if (requireNamespace("scater", quietly = TRUE)) {
          new_data <- scater::runUMAP(new_data, dimred = "harmony")
        }
      }

      # Store projected data
      projection_results$projected_data <- new_data
      projection_results$projected_embeddings <- projected_harmony

    } else if (integration_method == "mnn") {
      # For MNN integration
      if (!requireNamespace("batchelor", quietly = TRUE)) {
        stop("Package 'batchelor' required for projection with MNN integration")
      }

      # Project using batchelor
      projected_mnn <- batchelor::fastMNN(
        new_expr[common_genes, ],
        ref_expr[common_genes, ],
        d = 30
      )

      # Add projected dimensions
      SingleCellExperiment::reducedDim(new_data, "MNN") <- projected_mnn$corrected

      # Run UMAP on projected data
      if ("UMAP" %in% SingleCellExperiment::reducedDimNames(ref_data)) {
        if (requireNamespace("scater", quietly = TRUE)) {
          new_data <- scater::runUMAP(new_data, dimred = "MNN")
        }
      }

      # Store projected data
      projection_results$projected_data <- new_data
      projection_results$projected_embeddings <- projected_mnn$corrected
    }
  }

  # Create visualization of projection
  if (!is.null(projection_results$projected_data) && requireNamespace("ggplot2", quietly = TRUE)) {
    if (is_seurat) {
      # Extract embeddings
      if ("umap" %in% names(integrated_data@reductions)) {
        proj_embedding <- Seurat::Embeddings(projection_results$projected_data, reduction = "umap")
        orig_embedding <- Seurat::Embeddings(integrated_data, reduction = "umap")

        # Combine original and projected embeddings for visualization
        combined_df <- data.frame(
          rbind(orig_embedding, proj_embedding),
          Type = c(rep("Original", nrow(orig_embedding)), rep("Projected", nrow(proj_embedding))),
          Modality = c(metadata$modality, rep("New", nrow(proj_embedding)))
        )

        # Create plot
        dim_plot <- ggplot2::ggplot(combined_df, ggplot2::aes(x = UMAP_1, y = UMAP_2,
                                                              color = Modality, shape = Type)) +
          ggplot2::geom_point(size = 1, alpha = 0.7) +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = "Projection of New Data onto Integrated Space",
                        x = "UMAP1", y = "UMAP2")

        projection_results$plots$dim_plot <- dim_plot
      }
    } else {
      # Extract embeddings
      if ("UMAP" %in% SingleCellExperiment::reducedDimNames(integrated_data)) {
        proj_embedding <- SingleCellExperiment::reducedDim(projection_results$projected_data, "UMAP")
        orig_embedding <- SingleCellExperiment::reducedDim(integrated_data, "UMAP")

        # Combine original and projected embeddings for visualization
        combined_df <- data.frame(
          rbind(orig_embedding, proj_embedding),
          Type = c(rep("Original", nrow(orig_embedding)), rep("Projected", nrow(proj_embedding))),
          Modality = c(metadata$modality, rep("New", nrow(proj_embedding)))
        )

        # Create plot
        dim_plot <- ggplot2::ggplot(combined_df, ggplot2::aes(x = V1, y = V2,
                                                              color = Modality, shape = Type)) +
          ggplot2::geom_point(size = 1, alpha = 0.7) +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = "Projection of New Data onto Integrated Space",
                        x = "UMAP1", y = "UMAP2")

        projection_results$plots$dim_plot <- dim_plot
      }
    }
  }

  return(projection_results)
}

#' Helper function to project new data to existing low dimensional space
#'
#' @param new_data Matrix of new data in PCA space
#' @param ref_data Matrix of reference data in PCA space
#' @param ref_embedding Matrix of reference data in target low dimensional space
#'
#' @return Matrix of projected data in target low dimensional space
#' @keywords internal
projectToLowDim <- function(new_data, ref_data, ref_embedding) {
  # Find nearest neighbors in PCA space
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package 'FNN' required for projection")
  }

  # Ensure same number of dimensions
  n_dims <- min(ncol(new_data), ncol(ref_data))
  new_data <- new_data[, 1:n_dims, drop = FALSE]
  ref_data <- ref_data[, 1:n_dims, drop = FALSE]

  # Find k nearest neighbors for each new cell
  k <- min(20, nrow(ref_data))
  nn <- FNN::get.knnx(data = ref_data, query = new_data, k = k)

  # Calculate weights based on distances
  weights <- 1 / (nn$nn.dist + 1e-10)  # Add small value to avoid division by zero
  weights <- weights / rowSums(weights)  # Normalize weights

  # Project to low dimensional space as weighted average of neighbors
  projected <- matrix(0, nrow = nrow(new_data), ncol = ncol(ref_embedding))

  for (i in 1:nrow(new_data)) {
    neighbors <- nn$nn.index[i, ]
    projected[i, ] <- colSums(weights[i, ] * ref_embedding[neighbors, , drop = FALSE])
  }

  # Match rownames
  rownames(projected) <- rownames(new_data)
  colnames(projected) <- colnames(ref_embedding)

  return(projected)
}

#' Evaluate Integration Quality
#'
#' Evaluates the quality of multi-modal integration using various metrics.
#'
#' @param integration_results Output from integrateMultiModalData function
#' @param metrics Vector of metrics to compute (default: all available metrics)
#' @param ref_labels Optional reference cell labels for supervised evaluation
#' @param cell_type_var Name of the column containing cell type annotations (if available)
#' @param knn Number of nearest neighbors for graph-based metrics (default: 15)
#'
#' @return A list containing integration quality metrics and visualizations
#' @export
evaluateIntegration <- function(integration_results,
                                metrics = c("mixing", "silhouette", "kbet", "lisi"),
                                ref_labels = NULL,
                                cell_type_var = NULL,
                                knn = 15) {

  # Extract integrated data and embeddings
  integrated_data <- integration_results$integration_results
  embeddings <- integration_results$integration_embeddings

  # Check if input is valid
  if (is.null(integrated_data) || is.null(embeddings)) {
    stop("No integration results or embeddings found")
  }

  # Determine object type
  is_seurat <- inherits(integrated_data, "Seurat")
  is_sce <- inherits(integrated_data, "SingleCellExperiment")

  # Extract metadata
  if (is_seurat) {
    metadata <- integrated_data@meta.data
  } else {
    metadata <- as.data.frame(SingleCellExperiment::colData(integrated_data))
  }

  # Check if required columns exist in metadata
  if (!"modality" %in% colnames(metadata)) {
    stop("Column 'modality' not found in metadata")
  }

  # Initialize results
  eval_results <- list(
    metadata = list(
      metrics = metrics,
      knn = knn
    ),
    metrics = list(),
    plots = list()
  )

  # Extract modalities
  modalities <- unique(metadata$modality)

  # Check metric validity
  valid_metrics <- c("mixing", "silhouette", "kbet", "lisi", "ari", "nmi")
  invalid_metrics <- setdiff(metrics, valid_metrics)

  if (length(invalid_metrics) > 0) {
    warning(paste0("Invalid metrics: ", paste(invalid_metrics, collapse = ", "),
                   ". Will only compute valid metrics."))
    metrics <- intersect(metrics, valid_metrics)
  }

  # Compute KNN graph for graph-based metrics
  if (any(c("mixing", "kbet", "lisi") %in% metrics) && requireNamespace("FNN", quietly = TRUE)) {
    knn_results <- FNN::get.knn(embeddings, k = knn)

    # Create indices and distances matrices
    nn_idx <- knn_results$nn.index
    nn_dist <- knn_results$nn.dist

    # Set rownames
    rownames(nn_idx) <- rownames(embeddings)
    rownames(nn_dist) <- rownames(embeddings)
  }

  # Compute mixing metric
  if ("mixing" %in% metrics) {
    # Compute modality mixing score
    modality_mixing <- computeModalityMixing(nn_idx, metadata$modality)

    # Store results
    eval_results$metrics$mixing <- modality_mixing

    # Create visualization
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      # Compute summary statistics
      mixing_stats <- data.frame(
        Modality = names(modality_mixing$by_modality),
        MixingScore = unlist(modality_mixing$by_modality)
      )

      # Create barplot
      mixing_plot <- ggplot2::ggplot(mixing_stats, ggplot2::aes(x = Modality, y = MixingScore,
                                                                fill = Modality)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Modality Mixing Score",
                      subtitle = paste0("Overall: ", round(modality_mixing$overall, 3)),
                      x = "Modality", y = "Mixing Score")

      eval_results$plots$mixing <- mixing_plot
    }
  }

  # Compute silhouette metric
  if ("silhouette" %in% metrics && requireNamespace("cluster", quietly = TRUE)) {
    # Compute silhouette coefficient for modality separation
    # Convert modalities to numeric
    mod_ids <- as.integer(factor(metadata$modality))

    # Compute silhouette
    sil <- cluster::silhouette(mod_ids, stats::dist(embeddings))

    # Calculate summary statistics
    sil_avg <- mean(sil[, "sil_width"])
    sil_by_modality <- tapply(sil[, "sil_width"], metadata$modality, mean)

    # Store results
    eval_results$metrics$silhouette <- list(
      overall = sil_avg,
      by_modality = sil_by_modality,
      full = sil
    )

    # Create visualization
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      # Reshape data for plotting
      sil_df <- data.frame(
        Cell = rownames(embeddings),
        Modality = metadata$modality,
        Silhouette = sil[, "sil_width"]
      )

      # Create boxplot
      sil_plot <- ggplot2::ggplot(sil_df, ggplot2::aes(x = Modality, y = Silhouette, fill = Modality)) +
        ggplot2::geom_boxplot() +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Silhouette Coefficient by Modality",
                      subtitle = paste0("Overall: ", round(sil_avg, 3)),
                      x = "Modality", y = "Silhouette Width")

      eval_results$plots$silhouette <- sil_plot
    }
  }

  # Compute kBET metric (k-nearest neighbor batch effect test)
  if ("kbet" %in% metrics) {
    # Check if kBET package is available or implement simple version
    if (requireNamespace("kBET", quietly = TRUE)) {
      # Use kBET package
      kbet_results <- list()

      # Compute for each modality
      for (mod in modalities) {
        # Get cells from this modality
        mod_cells <- metadata$modality == mod

        # Skip if too few cells
        if (sum(mod_cells) < 10) {
          next
        }

        # Create batch vector (1 for this modality, 0 for others)
        batch <- as.factor(mod_cells)

        # Run kBET
        tryCatch({
          kbet_out <- kBET::kBET(embeddings, batch, k0 = min(knn, sum(mod_cells) - 1))
          kbet_results[[mod]] <- kbet_out$average.pval
        }, error = function(e) {
          warning(paste0("kBET failed for modality ", mod, ": ", e$message))
          kbet_results[[mod]] <- NA
        })
      }

      # Store results
      eval_results$metrics$kbet <- list(
        by_modality = kbet_results,
        overall = mean(unlist(kbet_results), na.rm = TRUE)
      )
    } else {
      # Implement simplified version of kBET
      kbet_simple <- computeSimpleKBET(nn_idx, metadata$modality)

      # Store results
      eval_results$metrics$kbet <- kbet_simple
    }

    # Create visualization
    if (requireNamespace("ggplot2", quietly = TRUE) && !is.null(eval_results$metrics$kbet)) {
      # Extract results
      if ("by_modality" %in% names(eval_results$metrics$kbet)) {
        kbet_stats <- data.frame(
          Modality = names(eval_results$metrics$kbet$by_modality),
          Acceptance = unlist(eval_results$metrics$kbet$by_modality)
        )
      } else {
        kbet_stats <- data.frame(
          Modality = names(eval_results$metrics$kbet$acceptance_rate),
          Acceptance = unlist(eval_results$metrics$kbet$acceptance_rate)
        )
      }

      # Create barplot
      kbet_plot <- ggplot2::ggplot(kbet_stats, ggplot2::aes(x = Modality, y = Acceptance,
                                                            fill = Modality)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme_minimal() +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(title = "kBET Acceptance Rate by Modality",
                      x = "Modality", y = "Acceptance Rate")

      eval_results$plots$kbet <- kbet_plot
    }
  }

  # Compute LISI metric (Local Inverse Simpson's Index)
  if ("lisi" %in% metrics) {
    # Check if LISI package is available or implement simple version
    if (requireNamespace("lisi", quietly = TRUE)) {
      # Use LISI package
      lisi_results <- lisi::compute_lisi(embeddings, metadata, "modality")

      # Store results
      eval_results$metrics$lisi <- list(
        scores = lisi_results,
        overall = mean(lisi_results$modality)
      )
    } else {
      # Implement simplified version of LISI
      lisi_simple <- computeSimpleLISI(nn_idx, metadata$modality)

      # Store results
      eval_results$metrics$lisi <- lisi_simple
    }

    # Create visualization
    if (requireNamespace("ggplot2", quietly = TRUE) && !is.null(eval_results$metrics$lisi)) {
      # Extract results
      if ("scores" %in% names(eval_results$metrics$lisi)) {
        lisi_df <- data.frame(
          Cell = rownames(embeddings),
          LISI = eval_results$metrics$lisi$scores$modality
        )
      } else {
        lisi_df <- data.frame(
          Cell = rownames(embeddings),
          LISI = eval_results$metrics$lisi$scores
        )
      }

      # Add modality information
      lisi_df$Modality <- metadata$modality

      # Create histogram
      lisi_plot <- ggplot2::ggplot(lisi_df, ggplot2::aes(x = LISI, fill = Modality)) +
        ggplot2::geom_histogram(position = "dodge", bins = 30) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Local Inverse Simpson's Index (LISI)",
                      subtitle = paste0("Mean: ", round(mean(lisi_df$LISI, na.rm = TRUE), 3)),
                      x = "LISI Score", y = "Count")

      eval_results$plots$lisi <- lisi_plot
    }
  }

  # Compute supervised metrics if reference labels are provided
  if (!is.null(ref_labels) || (!is.null(cell_type_var) && cell_type_var %in% colnames(metadata))) {
    # Get cell type labels
    if (is.null(ref_labels)) {
      cell_labels <- metadata[[cell_type_var]]
    } else {
      cell_labels <- ref_labels
    }

    # Compute Adjusted Rand Index (ARI)
    if ("ari" %in% metrics && requireNamespace("mclust", quietly = TRUE)) {
      # Cluster integrated data
      if (requireNamespace("stats", quietly = TRUE)) {
        # Hierarchical clustering
        hc <- stats::hclust(stats::dist(embeddings), method = "ward.D2")
        clusters <- stats::cutree(hc, k = length(unique(cell_labels)))

        # Compute ARI
        ari <- mclust::adjustedRandIndex(cell_labels, clusters)

        # Store results
        eval_results$metrics$ari <- ari
      }
    }

    # Compute Normalized Mutual Information (NMI)
    if ("nmi" %in% metrics && requireNamespace("aricode", quietly = TRUE)) {
      # Cluster integrated data if not already done
      if (!exists("clusters")) {
        # Hierarchical clustering
        hc <- stats::hclust(stats::dist(embeddings), method = "ward.D2")
        clusters <- stats::cutree(hc, k = length(unique(cell_labels)))
      }

      # Compute NMI
      nmi <- aricode::NMI(cell_labels, clusters)

      # Store results
      eval_results$metrics$nmi <- nmi
    }

    # Create summary plot for supervised metrics
    if (requireNamespace("ggplot2", quietly = TRUE) &&
        (("ari" %in% names(eval_results$metrics)) || ("nmi" %in% names(eval_results$metrics)))) {

      # Create data frame for plotting
      metrics_df <- data.frame(
        Metric = character(),
        Value = numeric(),
        stringsAsFactors = FALSE
      )

      if ("ari" %in% names(eval_results$metrics)) {
        metrics_df <- rbind(metrics_df, data.frame(
          Metric = "ARI",
          Value = eval_results$metrics$ari
        ))
      }

      if ("nmi" %in% names(eval_results$metrics)) {
        metrics_df <- rbind(metrics_df, data.frame(
          Metric = "NMI",
          Value = eval_results$metrics$nmi
        ))
      }

      # Create barplot
      sup_plot <- ggplot2::ggplot(metrics_df, ggplot2::aes(x = Metric, y = Value, fill = Metric)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme_minimal() +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(title = "Supervised Evaluation Metrics",
                      x = "", y = "Score")

      eval_results$plots$supervised <- sup_plot
    }
  }

  # Return results
  return(eval_results)
}

#' Compute modality mixing score
#'
#' @param nn_idx KNN indices matrix
#' @param modality_labels Vector of modality labels
#'
#' @return List with mixing scores
#' @keywords internal
computeModalityMixing <- function(nn_idx, modality_labels) {
  # Get unique modalities
  modalities <- unique(modality_labels)

  # Initialize mixing scores
  mixing_scores <- numeric(nrow(nn_idx))
  names(mixing_scores) <- rownames(nn_idx)

  # Calculate mixing score for each cell
  for (i in 1:nrow(nn_idx)) {
    # Get modalities of neighbors
    neighbor_mods <- modality_labels[nn_idx[i, ]]

    # Calculate mixing score (proportion of neighbors from different modality)
    cell_mod <- modality_labels[i]
    mixing_scores[i] <- mean(neighbor_mods != cell_mod)
  }

  # Calculate average mixing score by modality
  mixing_by_modality <- tapply(mixing_scores, modality_labels, mean)

  # Return results
  return(list(
    overall = mean(mixing_scores),
    by_modality = mixing_by_modality,
    scores = mixing_scores
  ))
}

#' Compute simplified kBET metric
#'
#' @param nn_idx KNN indices matrix
#' @param modality_labels Vector of modality labels
#'
#' @return List with kBET results
#' @keywords internal
computeSimpleKBET <- function(nn_idx, modality_labels) {
  # Get unique modalities
  modalities <- unique(modality_labels)
  n_modalities <- length(modalities)

  # Get expected frequency of each modality
  expected_freq <- table(modality_labels) / length(modality_labels)

  # Initialize results
  kbet_scores <- numeric(nrow(nn_idx))
  names(kbet_scores) <- rownames(nn_idx)

  # Calculate kBET score for each cell
  for (i in 1:nrow(nn_idx)) {
    # Get modalities of neighbors
    neighbor_mods <- modality_labels[nn_idx[i, ]]

    # Calculate observed frequency
    observed_freq <- table(factor(neighbor_mods, levels = modalities)) / length(neighbor_mods)

    # Calculate chi-square statistic
    chi_sq <- sum(((observed_freq - expected_freq)^2) / expected_freq)

    # Convert to p-value
    p_val <- 1 - stats::pchisq(chi_sq, df = n_modalities - 1)

    # Store p-value as kBET score
    kbet_scores[i] <- p_val
  }

  # Calculate acceptance rate (proportion of cells with p-value > 0.05)
  acceptance_rate <- tapply(kbet_scores, modality_labels, function(x) mean(x > 0.05))

  # Return results
  return(list(
    scores = kbet_scores,
    acceptance_rate = acceptance_rate,
    overall = mean(kbet_scores > 0.05)
  ))
}

#' Compute simplified LISI metric
#'
#' @param nn_idx KNN indices matrix
#' @param modality_labels Vector of modality labels
#'
#' @return List with LISI results
#' @keywords internal
computeSimpleLISI <- function(nn_idx, modality_labels) {
  # Initialize results
  lisi_scores <- numeric(nrow(nn_idx))
  names(lisi_scores) <- rownames(nn_idx)

  # Calculate LISI score for each cell
  for (i in 1:nrow(nn_idx)) {
    # Get modalities of neighbors
    neighbor_mods <- modality_labels[nn_idx[i, ]]

    # Calculate Simpson index (probability that two randomly selected neighbors are from the same modality)
    mod_freq <- table(neighbor_mods) / length(neighbor_mods)
    simpson_index <- sum(mod_freq^2)

    # Inverse Simpson index (effective number of modalities in the neighborhood)
    lisi_scores[i] <- 1 / simpson_index
  }

  # Calculate statistics by modality
  lisi_by_modality <- tapply(lisi_scores, modality_labels, mean)

  # Return results
  return(list(
    scores = lisi_scores,
    by_modality = lisi_by_modality,
    overall = mean(lisi_scores)
  ))
}
#' Apply Integrated Analysis to Identify Cell Type-Specific Multi-Modal Signatures
#'
#' Identifies cell type-specific multi-modal signatures by integrating information
#' across different data modalities.
#'
#' @param integration_results Output from integrateMultiModalData function
#' @param cell_type_var Name of the column containing cell type annotations
#' @param n_features Number of top features to include in signatures (default: 50)
#' @param min_cells Minimum number of cells required per cell type (default: 10)
#' @param concordance_threshold Minimum concordance score to include features (default: 0.7)
#' @param modality_weights Optional named vector of weights for each modality
#' @param return_models Logical, whether to return trained classification models (default: FALSE)
#'
#' @return A list containing multi-modal signatures for each cell type
#' @export
findMultiModalSignatures <- function(integration_results,
                                     cell_type_var,
                                     n_features = 50,
                                     min_cells = 10,
                                     concordance_threshold = 0.7,
                                     modality_weights = NULL,
                                     return_models = FALSE) {

  # Extract integrated data
  integrated_data <- integration_results$integration_results

  # Check if input is valid
  if (is.null(integrated_data)) {
    stop("No integration results found")
  }

  # Determine object type
  is_seurat <- inherits(integrated_data, "Seurat")
  is_sce <- inherits(integrated_data, "SingleCellExperiment")

  if (!is_seurat && !is_sce) {
    stop("Integrated data must be a Seurat or SingleCellExperiment object")
  }

  # Extract metadata
  if (is_seurat) {
    metadata <- integrated_data@meta.data
  } else {
    metadata <- as.data.frame(SingleCellExperiment::colData(integrated_data))
  }

  # Check if required columns exist in metadata
  if (!cell_type_var %in% colnames(metadata)) {
    stop(paste0("Column '", cell_type_var, "' not found in metadata"))
  }

  if (!"modality" %in% colnames(metadata)) {
    stop("Column 'modality' not found in metadata")
  }

  # Extract cell types and modalities
  cell_types <- unique(metadata[[cell_type_var]])
  modalities <- unique(metadata$modality)

  # Validate modality weights if provided
  if (!is.null(modality_weights)) {
    if (!all(names(modality_weights) %in% modalities)) {
      stop("Names of modality_weights must match modalities in the data")
    }
  } else {
    # Set equal weights
    modality_weights <- rep(1, length(modalities))
    names(modality_weights) <- modalities
  }

  # Initialize results
  signature_results <- list(
    metadata = list(
      cell_types = cell_types,
      modalities = modalities,
      modality_weights = modality_weights,
      n_features = n_features,
      concordance_threshold = concordance_threshold
    ),
    signatures = list(),
    feature_stats = list(),
    plots = list()
  )

  # For each cell type, identify multi-modal signature
  for (cell_type in cell_types) {
    # Get cells of this cell type
    cell_type_cells <- rownames(metadata)[metadata[[cell_type_var]] == cell_type]

    # Skip if too few cells
    if (length(cell_type_cells) < min_cells) {
      message(paste0("Skipping cell type ", cell_type, ": insufficient cells (",
                     length(cell_type_cells), " < ", min_cells, ")"))
      next
    }

    # Identify differentially expressed features in each modality
    modality_markers <- list()
    modality_models <- list()

    for (modality in modalities) {
      # Get cells from this modality
      modality_cells <- rownames(metadata)[metadata$modality == modality]

      # Skip if too few cells
      if (length(modality_cells) < min_cells) {
        message(paste0("Skipping modality ", modality, " for cell type ", cell_type,
                       ": insufficient cells"))
        next
      }

      # Subset data to this modality
      if (is_seurat) {
        modality_data <- subset(integrated_data, cells = modality_cells)

        # Set cell type identity
        modality_data$cell_type_binary <- modality_data@meta.data[[cell_type_var]] == cell_type
        Seurat::Idents(modality_data) <- "cell_type_binary"

        # Find markers
        markers <- Seurat::FindMarkers(
          modality_data,
          ident.1 = TRUE,  # This cell type
          ident.2 = FALSE, # Other cell types
          logfc.threshold = 0.25,
          test.use = "wilcox",
          min.pct = 0.1
        )

        # Add gene names
        markers$feature <- rownames(markers)

      } else {
        # For SingleCellExperiment
        modality_data <- integrated_data[, modality_cells]

        # Create cell type binary vector
        cell_type_binary <- SingleCellExperiment::colData(modality_data)[[cell_type_var]] == cell_type

        # Find markers using scran
        if (requireNamespace("scran", quietly = TRUE)) {
          design <- model.matrix(~cell_type_binary)
          markers_result <- scran::findMarkers(
            modality_data,
            groups = cell_type_binary,
            design = design,
            direction = "up"
          )

          # Format results
          markers <- as.data.frame(markers_result[[TRUE]])
          markers$feature <- rownames(markers)
          markers$p_val <- markers$p.value
          markers$p_val_adj <- markers$FDR
          markers$avg_log2FC <- markers$summary.logFC
        } else {
          # Simple implementation if scran not available
          if ("logcounts" %in% SingleCellExperiment::assayNames(modality_data)) {
            expr <- SingleCellExperiment::logcounts(modality_data)
          } else {
            expr <- SingleCellExperiment::counts(modality_data)
          }

          # Calculate log fold change
          avg_expr_in <- rowMeans(expr[, cell_type_binary])
          avg_expr_out <- rowMeans(expr[, !cell_type_binary])
          logfc <- avg_expr_in - avg_expr_out

          # Calculate p-values using t-test
          p_vals <- sapply(1:nrow(expr), function(i) {
            t_result <- tryCatch({
              stats::t.test(expr[i, cell_type_binary], expr[i, !cell_type_binary])$p.value
            }, error = function(e) {
              1  # Return 1 if t-test fails
            })
            return(t_result)
          })

          # Adjust p-values
          p_adj <- stats::p.adjust(p_vals, method = "BH")

          # Combine results
          markers <- data.frame(
            feature = rownames(expr),
            avg_log2FC = logfc,
            p_val = p_vals,
            p_val_adj = p_adj,
            stringsAsFactors = FALSE
          )
        }
      }

      # Filter markers and add modality information
      markers$modality <- modality
      markers$significance <- -log10(markers$p_val_adj)
      markers$direction <- ifelse(markers$avg_log2FC > 0, "up", "down")

      # Keep only upregulated markers with adjusted p-value < 0.05
      filtered_markers <- markers[markers$p_val_adj < 0.05 & markers$avg_log2FC > 0, ]

      # Store results
      modality_markers[[modality]] <- filtered_markers

      # Train classifier model if requested
      if (return_models) {
        # Extract expression data
        if (is_seurat) {
          expr_matrix <- t(as.matrix(Seurat::GetAssayData(modality_data, slot = "data")))
        } else {
          if ("logcounts" %in% SingleCellExperiment::assayNames(modality_data)) {
            expr_matrix <- t(as.matrix(SingleCellExperiment::logcounts(modality_data)))
          } else {
            expr_matrix <- t(as.matrix(SingleCellExperiment::counts(modality_data)))
          }
        }

        # Create training data
        training_data <- data.frame(
          cell_type = cell_type_binary,
          expr_matrix,
          stringsAsFactors = FALSE
        )

        # Train model using randomForest if available
        if (requireNamespace("randomForest", quietly = TRUE)) {
          # Select top features
          top_marker_features <- filtered_markers$feature[order(filtered_markers$p_val_adj)][1:min(50, nrow(filtered_markers))]
          training_features <- intersect(top_marker_features, colnames(training_data))

          if (length(training_features) > 0) {
            model_formula <- stats::as.formula(paste("cell_type ~",
                                                     paste(training_features, collapse = " + ")))

            rf_model <- tryCatch({
              randomForest::randomForest(model_formula, data = training_data, ntree = 100)
            }, error = function(e) {
              NULL
            })

            modality_models[[modality]] <- rf_model
          }
        }
      }
    }

    # Integrate markers across modalities
    all_markers <- do.call(rbind, modality_markers)

    if (!is.null(all_markers) && nrow(all_markers) > 0) {
      # Weight markers by modality
      weighted_markers <- all_markers
      for (modality in modalities) {
        if (modality %in% names(modality_weights)) {
          weighted_markers$significance[weighted_markers$modality == modality] <-
            weighted_markers$significance[weighted_markers$modality == modality] * modality_weights[modality]
        }
      }

      # Calculate concordance score for each feature
      # Features appearing in multiple modalities will have higher scores
      feature_counts <- table(weighted_markers$feature)
      modality_counts <- table(weighted_markers$modality)

      # Map feature count to each marker
      weighted_markers$n_modalities <- feature_counts[weighted_markers$feature]

      # Calculate concordance score (proportion of modalities where the feature is a marker)
      weighted_markers$concordance <- weighted_markers$n_modalities / length(modalities)

      # Filter for high concordance markers
      high_concordance <- weighted_markers[weighted_markers$concordance >= concordance_threshold, ]

      # Collapse by feature - keep the highest significance score
      if (nrow(high_concordance) > 0) {
        feature_summary <- stats::aggregate(
          significance ~ feature,
          data = high_concordance,
          FUN = max
        )

        # Add concordance information
        feature_summary$concordance <- feature_counts[feature_summary$feature] / length(modalities)

        # Sort by significance and concordance
        feature_summary <- feature_summary[order(-feature_summary$concordance, -feature_summary$significance), ]

        # Get top features as signature
        signature_features <- feature_summary$feature[1:min(n_features, nrow(feature_summary))]

        # Store signature
        signature_results$signatures[[cell_type]] <- signature_features
        signature_results$feature_stats[[cell_type]] <- feature_summary

        # Create visualization if ggplot2 is available
        if (requireNamespace("ggplot2", quietly = TRUE) && nrow(feature_summary) > 0) {
          # Plot top features
          top_n_for_plot <- min(20, nrow(feature_summary))
          plot_data <- feature_summary[1:top_n_for_plot, ]

          sig_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = reorder(feature, significance),
                                                              y = significance,
                                                              fill = concordance)) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::scale_fill_viridis_c(name = "Concordance") +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = paste0("Multi-Modal Signature for ", cell_type),
                          subtitle = paste0(length(signature_features), " features"),
                          x = "Feature", y = "-log10(p-value)") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

          signature_results$plots[[cell_type]] <- sig_plot
        }
      }
    }

    # Add models if requested
    if (return_models) {
      signature_results$models[[cell_type]] <- modality_models
    }
  }

  # Create summary visualization
  if (requireNamespace("ggplot2", quietly = TRUE) && length(signature_results$signatures) > 0) {
    # Create summary data
    summary_data <- data.frame(
      CellType = names(signature_results$signatures),
      SignatureSize = sapply(signature_results$signatures, length)
    )

    # Plot summary
    summary_plot <- ggplot2::ggplot(summary_data, ggplot2::aes(x = reorder(CellType, SignatureSize),
                                                               y = SignatureSize,
                                                               fill = CellType)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Multi-Modal Signature Summary",
                    x = "Cell Type", y = "Number of Features") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     legend.position = "none")

    signature_results$plots$summary <- summary_plot
  }

  return(signature_results)
}
