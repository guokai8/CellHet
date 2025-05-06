#' Optimize Parameters for Multi-Modal Integration
#'
#' Automatically tests different integration parameters and selects optimal settings
#' based on quality metrics, helping users achieve the best possible integration results.
#'
#' @param data_list List of data objects (Seurat or SingleCellExperiment) for different modalities
#' @param modality_names Names of modalities corresponding to each item in data_list
#' @param cell_type_var Name of the variable containing cell type annotations
#' @param optimization_metric Primary metric to optimize: "mixing", "silhouette", "kbet", or "lisi" (default: "mixing")
#' @param secondary_metrics Additional metrics to compute (default: all available metrics)
#' @param methods Vector of integration methods to test (default: all available methods)
#' @param n_dims_range Range of dimensions to test (default: c(10, 20, 30))
#' @param n_features_range Range of features to test (default: c(1000, 2000, 3000))
#' @param max_combinations Maximum number of parameter combinations to test (default: 10)
#' @param n_iterations Number of iterations for more reliable metrics (default: 1)
#' @param reference_batch Optional reference batch/modality for comparison
#' @param return_all_results Logical, whether to return all results or just the best (default: FALSE)
#' @param parallel Logical, whether to use parallel processing (default: FALSE)
#' @param n_cores Number of cores to use for parallel processing (default: 2)
#' @param seed Random seed for reproducibility (default: 42)
#' @param verbose Logical, whether to show progress messages (default: TRUE)
#'
#' @return A list containing optimized parameters and integration results
#' @export
optimizeIntegrationParameters <- function(data_list,
                                          modality_names = NULL,
                                          cell_type_var = "cell_type",
                                          optimization_metric = "mixing",
                                          secondary_metrics = c("silhouette", "kbet", "lisi"),
                                          methods = c("canonical_correlation", "harmony", "mnn", "seurat"),
                                          n_dims_range = c(10, 20, 30),
                                          n_features_range = c(1000, 2000, 3000),
                                          max_combinations = 10,
                                          n_iterations = 1,
                                          reference_batch = NULL,
                                          return_all_results = FALSE,
                                          parallel = FALSE,
                                          n_cores = 2,
                                          seed = 42,
                                          verbose = TRUE) {

  # Set seed for reproducibility
  set.seed(seed)

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

  # Validate that all objects are of the same type
  if (length(unique(object_types)) > 1) {
    stop("All objects must be of the same type (either all Seurat or all SingleCellExperiment)")
  }

  # Validate optimization metric
  valid_metrics <- c("mixing", "silhouette", "kbet", "lisi")
  if (!optimization_metric %in% valid_metrics) {
    stop(paste0("Invalid optimization metric. Must be one of: ", paste(valid_metrics, collapse = ", ")))
  }

  # Filter methods based on object type
  if (unique(object_types) == "SCE") {
    methods <- methods[methods %in% c("canonical_correlation", "harmony", "mnn")]
    if (length(methods) == 0) {
      stop("No valid integration methods for SingleCellExperiment objects")
    }
    if ("seurat" %in% methods) {
      warning("Seurat integration method not available for SingleCellExperiment objects. Removing from methods.")
      methods <- methods[methods != "seurat"]
    }
  }

  # Generate parameter grid
  param_grid <- expand.grid(
    method = methods,
    n_dims = n_dims_range,
    n_features = n_features_range,
    stringsAsFactors = FALSE
  )

  # Limit number of combinations if needed
  if (nrow(param_grid) > max_combinations) {
    if (verbose) {
      message(paste0("Limiting parameter combinations from ", nrow(param_grid), " to ", max_combinations))
    }
    param_grid <- param_grid[sample(1:nrow(param_grid), max_combinations), ]
  }

  if (verbose) {
    message(paste0("Testing ", nrow(param_grid), " parameter combinations"))
  }

  # Set up parallelization if requested
  if (parallel) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("Package 'parallel' required for parallel processing. Falling back to sequential.")
      parallel <- FALSE
    } else {
      if (verbose) {
        message(paste0("Setting up parallel processing with ", n_cores, " cores"))
      }
      cl <- parallel::makeCluster(n_cores)
      parallel::clusterExport(cl, c("data_list", "modality_names", "cell_type_var",
                                    "optimization_metric", "secondary_metrics",
                                    "n_iterations", "reference_batch", "seed", "verbose"))

      # Load necessary packages on all workers
      parallel::clusterEvalQ(cl, {
        library(methods)
        if (requireNamespace("Seurat", quietly = TRUE)) library(Seurat)
        if (requireNamespace("SingleCellExperiment", quietly = TRUE)) library(SingleCellExperiment)
      })
    }
  }

  # Function to evaluate a single parameter combination
  evaluateParams <- function(params_row) {
    # Extract parameters
    method <- params_row["method"]
    n_dims <- as.numeric(params_row["n_dims"])
    n_features <- as.numeric(params_row["n_features"])

    # Initialize results for this parameter set
    param_results <- list(
      parameters = list(
        method = method,
        n_dims = n_dims,
        n_features = n_features
      ),
      integration_results = NULL,
      metrics = list()
    )

    # Perform integration
    tryCatch({
      if (verbose) {
        message(paste0("Testing: method=", method, ", n_dims=", n_dims, ", n_features=", n_features))
      }

      # Run integration
      integration_result <- integrateMultiModalData(
        data_list = data_list,
        modality_names = modality_names,
        cell_type_var = cell_type_var,
        integration_method = method,
        dim_reduction = "umap",
        n_dims = n_dims,
        n_features = n_features,
        seed = seed
      )

      # Store integration result
      param_results$integration_results <- integration_result

      # Extract integrated data and embeddings
      integrated_data <- integration_result$integration_results
      integrated_embeddings <- integration_result$integration_embeddings

      # Extract metadata
      if (inherits(integrated_data, "Seurat")) {
        metadata <- integrated_data@meta.data
      } else {
        metadata <- as.data.frame(SingleCellExperiment::colData(integrated_data))
      }

      # Compute integration quality metrics
      # 1. Modality mixing score
      if ("mixing" %in% c(optimization_metric, secondary_metrics)) {
        if (requireNamespace("FNN", quietly = TRUE)) {
          knn_results <- FNN::get.knn(integrated_embeddings, k = 15)
          nn_idx <- knn_results$nn.index

          # Compute modality mixing
          modality_mixing <- computeModalityMixing(nn_idx, metadata$modality)
          param_results$metrics$mixing <- modality_mixing$overall
        }
      }

      # 2. Silhouette score
      if ("silhouette" %in% c(optimization_metric, secondary_metrics)) {
        if (requireNamespace("cluster", quietly = TRUE)) {
          # Lower silhouette is better for batch integration (less separation between batches)
          mod_ids <- as.integer(factor(metadata$modality))
          sil <- cluster::silhouette(mod_ids, stats::dist(integrated_embeddings))
          param_results$metrics$silhouette <- mean(sil[, "sil_width"])
        }
      }

      # 3. kBET metric
      if ("kbet" %in% c(optimization_metric, secondary_metrics)) {
        if (requireNamespace("FNN", quietly = TRUE)) {
          knn_results <- FNN::get.knn(integrated_embeddings, k = 15)
          nn_idx <- knn_results$nn.index

          # Compute kBET
          kbet_simple <- computeSimpleKBET(nn_idx, metadata$modality)
          param_results$metrics$kbet <- kbet_simple$overall
        }
      }

      # 4. LISI metric
      if ("lisi" %in% c(optimization_metric, secondary_metrics)) {
        if (requireNamespace("FNN", quietly = TRUE)) {
          knn_results <- FNN::get.knn(integrated_embeddings, k = 15)
          nn_idx <- knn_results$nn.index

          # Compute LISI
          lisi_simple <- computeSimpleLISI(nn_idx, metadata$modality)
          param_results$metrics$lisi <- lisi_simple$overall
        }
      }

      return(param_results)
    }, error = function(e) {
      if (verbose) {
        message(paste0("Error with parameters: method=", method, ", n_dims=", n_dims, ", n_features=", n_features))
        message(paste0("Error message: ", e$message))
      }

      # Return partial results
      param_results$error <- e$message
      return(param_results)
    })
  }

  # Evaluate all parameter combinations
  if (parallel) {
    all_results <- parallel::parLapply(cl, split(param_grid, 1:nrow(param_grid)), evaluateParams)
    parallel::stopCluster(cl)
  } else {
    all_results <- lapply(split(param_grid, 1:nrow(param_grid)), evaluateParams)
  }

  # Extract metrics for ranking
  metrics_df <- data.frame(
    method = sapply(all_results, function(x) x$parameters$method),
    n_dims = sapply(all_results, function(x) x$parameters$n_dims),
    n_features = sapply(all_results, function(x) x$parameters$n_features),
    error = sapply(all_results, function(x) ifelse(is.null(x$error), FALSE, TRUE))
  )

  # Add metrics
  for (metric in c(optimization_metric, secondary_metrics)) {
    metrics_df[[metric]] <- sapply(all_results, function(x) {
      if (!is.null(x$metrics[[metric]])) {
        return(x$metrics[[metric]])
      } else {
        return(NA)
      }
    })
  }

  # Remove rows with errors
  metrics_df <- metrics_df[!metrics_df$error, ]

  # Rank parameter combinations based on optimization metric
  if (nrow(metrics_df) > 0) {
    # Determine ranking direction based on metric
    if (optimization_metric %in% c("mixing", "kbet")) {
      # Higher is better
      ranked_metrics <- metrics_df[order(-metrics_df[[optimization_metric]]), ]
    } else if (optimization_metric == "silhouette") {
      # Lower is better (less separation between batches)
      ranked_metrics <- metrics_df[order(metrics_df[[optimization_metric]]), ]
    } else if (optimization_metric == "lisi") {
      # Closer to number of batches is better
      n_batches <- length(unique(modality_names))
      ranked_metrics <- metrics_df[order(abs(metrics_df[[optimization_metric]] - n_batches)), ]
    }

    # Get best parameters
    best_params <- ranked_metrics[1, c("method", "n_dims", "n_features")]

    # Rerun integration with best parameters for final result
    if (verbose) {
      message(paste0("Best parameters: method=", best_params$method,
                     ", n_dims=", best_params$n_dims,
                     ", n_features=", best_params$n_features))
    }

    best_integration <- integrateMultiModalData(
      data_list = data_list,
      modality_names = modality_names,
      cell_type_var = cell_type_var,
      integration_method = best_params$method,
      dim_reduction = "umap",
      n_dims = best_params$n_dims,
      n_features = best_params$n_features,
      seed = seed
    )

    # Compile results
    optimization_results <- list(
      best_parameters = best_params,
      best_integration = best_integration,
      all_metrics = ranked_metrics
    )

    # Add all results if requested
    if (return_all_results) {
      optimization_results$all_results <- all_results
    }

    return(optimization_results)

  } else {
    stop("No successful parameter combinations found. Try different parameter ranges.")
  }
}

#' Integration Feature Importance Analysis
#'
#' Identifies and ranks features by their importance in distinguishing cell types
#' across modalities using machine learning approaches.
#'
#' @param integration_results Output from integrateMultiModalData function
#' @param cell_type_var Name of the column containing cell type annotations
#' @param n_features Number of top features to return (default: 50)
#' @param method Method for feature importance: "random_forest", "lasso", or "elastic_net" (default: "random_forest")
#' @param modality_specific Logical, whether to calculate feature importance for each modality separately (default: TRUE)
#' @param test_fraction Fraction of data to use for testing (default: 0.3)
#' @param n_iterations Number of iterations for stability analysis (default: 5)
#' @param balance_classes Logical, whether to balance classes for training (default: TRUE)
#' @param seed Random seed for reproducibility (default: 42)
#'
#' @return A list containing feature importance scores and visualizations
#' @export
analyzeFeatureImportance <- function(integration_results,
                                     cell_type_var,
                                     n_features = 50,
                                     method = "random_forest",
                                     modality_specific = TRUE,
                                     test_fraction = 0.3,
                                     n_iterations = 5,
                                     balance_classes = TRUE,
                                     seed = 42) {

  # Set seed for reproducibility
  set.seed(seed)

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
    warning("Column 'modality' not found in metadata. Setting modality_specific to FALSE.")
    modality_specific <- FALSE
  }

  # Extract expression data
  if (is_seurat) {
    # Get data from default assay
    expr_data <- as.matrix(Seurat::GetAssayData(integrated_data, slot = "data"))
  } else {
    # For SCE
    if ("logcounts" %in% SummarizedExperiment::assayNames(integrated_data)) {
      expr_data <- SummarizedExperiment::assay(integrated_data, "logcounts")
    } else {
      expr_data <- SummarizedExperiment::assay(integrated_data, 1)
    }
  }

  # Initialize results
  importance_results <- list(
    metadata = list(
      method = method,
      n_features = n_features,
      modality_specific = modality_specific,
      n_iterations = n_iterations
    ),
    feature_importance = list(),
    top_features = list(),
    stability = list(),
    performance = list(),
    plots = list()
  )

  # Function to run feature importance analysis
  runImportanceAnalysis <- function(expr_matrix, cell_labels, method) {
    # Remove genes with zero variance
    gene_vars <- apply(expr_matrix, 1, stats::var)
    valid_genes <- names(gene_vars[gene_vars > 0])
    expr_matrix <- expr_matrix[valid_genes, ]

    # Transpose for ML framework (samples x features)
    expr_matrix <- t(expr_matrix)

    # Split into training and test sets
    n_samples <- nrow(expr_matrix)
    test_indices <- sample(1:n_samples, size = round(n_samples * test_fraction))
    train_indices <- setdiff(1:n_samples, test_indices)

    # Create training and test sets
    train_data <- expr_matrix[train_indices, ]
    test_data <- expr_matrix[test_indices, ]
    train_labels <- cell_labels[train_indices]
    test_labels <- cell_labels[test_indices]

    # Balance classes if requested
    if (balance_classes) {
      # Get class frequencies
      class_counts <- table(train_labels)
      min_class_size <- min(class_counts)

      # Create balanced set
      balanced_indices <- c()
      for (class_name in names(class_counts)) {
        class_indices <- which(train_labels == class_name)
        if (length(class_indices) > min_class_size) {
          class_indices <- sample(class_indices, min_class_size)
        }
        balanced_indices <- c(balanced_indices, class_indices)
      }

      # Update training data
      train_data <- train_data[balanced_indices, ]
      train_labels <- train_labels[balanced_indices]
    }

    # Calculate feature importance based on method
    if (method == "random_forest") {
      if (requireNamespace("randomForest", quietly = TRUE)) {
        # Train random forest
        rf_model <- randomForest::randomForest(
          x = train_data,
          y = factor(train_labels),
          importance = TRUE,
          ntree = 100
        )

        # Extract feature importance
        importance_values <- randomForest::importance(rf_model)

        # For multi-class, use MeanDecreaseGini
        if (ncol(importance_values) > 2) {
          feature_importance <- importance_values[, "MeanDecreaseGini"]
        } else {
          feature_importance <- importance_values[, "MeanDecreaseAccuracy"]
        }

        # Make predictions on test set
        predictions <- predict(rf_model, test_data)

        # Calculate accuracy
        accuracy <- sum(predictions == test_labels) / length(test_labels)

        # Return results
        return(list(
          feature_importance = feature_importance,
          model = rf_model,
          performance = list(accuracy = accuracy)
        ))
      } else {
        stop("Package 'randomForest' required for random_forest method")
      }
    } else if (method == "lasso" || method == "elastic_net") {
      if (requireNamespace("glmnet", quietly = TRUE)) {
        # Set alpha based on method
        alpha_val <- ifelse(method == "lasso", 1, 0.5)  # 1 for lasso, 0.5 for elastic net

        # Create model matrix
        x_train <- train_data
        y_train <- train_labels

        # Train model
        cv_fit <- glmnet::cv.glmnet(
          x = x_train,
          y = factor(y_train),
          family = "multinomial",
          alpha = alpha_val,
          nfolds = 5
        )

        # Extract feature importance
        coef_matrix <- glmnet::coef.glmnet(cv_fit, s = "lambda.min")

        # Combine coefficients across classes
        feature_importance <- numeric(ncol(train_data))
        names(feature_importance) <- colnames(train_data)

        for (i in 1:length(coef_matrix)) {
          class_coef <- as.matrix(coef_matrix[[i]])
          feature_importance <- feature_importance + abs(class_coef[-1, 1])  # Skip intercept
        }

        # Make predictions on test set
        x_test <- test_data
        predictions <- predict(cv_fit, newx = x_test, s = "lambda.min", type = "class")

        # Calculate accuracy
        accuracy <- sum(predictions == test_labels) / length(test_labels)

        # Return results
        return(list(
          feature_importance = feature_importance,
          model = cv_fit,
          performance = list(accuracy = accuracy)
        ))
      } else {
        stop("Package 'glmnet' required for lasso or elastic_net method")
      }
    } else {
      stop(paste0("Unsupported method: ", method))
    }
  }

  # Run feature importance analysis
  if (modality_specific) {
    # Run for each modality separately
    modalities <- unique(metadata$modality)

    for (modality in modalities) {
      if (verbose) {
        message(paste0("Running feature importance analysis for modality: ", modality))
      }

      # Get cells from this modality
      modality_cells <- rownames(metadata)[metadata$modality == modality]

      # Extract cell types
      cell_types <- metadata[modality_cells, cell_type_var]

      # Extract expression data
      modality_expr <- expr_data[, modality_cells]

      # Run multiple iterations
      iteration_results <- list()
      for (i in 1:n_iterations) {
        iter_result <- runImportanceAnalysis(modality_expr, cell_types, method)
        iteration_results[[i]] <- iter_result
      }

      # Combine results across iterations
      combined_importance <- Reduce(`+`, lapply(iteration_results, function(x) x$feature_importance))
      combined_importance <- combined_importance / n_iterations

      # Get top features
      top_features_idx <- order(combined_importance, decreasing = TRUE)[1:min(n_features, length(combined_importance))]
      top_features <- names(combined_importance)[top_features_idx]

      # Calculate stability
      stability_scores <- numeric(length(top_features))
      names(stability_scores) <- top_features

      for (feature in top_features) {
        # Count how many times feature appears in top features across iterations
        feature_ranks <- sapply(iteration_results, function(x) {
          imp <- x$feature_importance
          which(names(imp)[order(imp, decreasing = TRUE)] == feature)
        })

        # Calculate stability as inverse of rank variance
        stability_scores[feature] <- 1 / (1 + stats::var(feature_ranks))
      }

      # Calculate average performance
      avg_accuracy <- mean(sapply(iteration_results, function(x) x$performance$accuracy))

      # Store results
      importance_results$feature_importance[[modality]] <- combined_importance
      importance_results$top_features[[modality]] <- top_features
      importance_results$stability[[modality]] <- stability_scores
      importance_results$performance[[modality]] <- list(accuracy = avg_accuracy)

      # Create visualization
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        # Create data frame for plotting
        plot_data <- data.frame(
          Feature = names(combined_importance[top_features_idx]),
          Importance = combined_importance[top_features_idx],
          Stability = stability_scores[names(combined_importance[top_features_idx])],
          stringsAsFactors = FALSE
        )

        # Create importance plot
        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = reorder(Feature, Importance), y = Importance, fill = Stability)) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::scale_fill_viridis_c() +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
          ggplot2::labs(
            title = paste0("Feature Importance for ", modality),
            subtitle = paste0("Prediction Accuracy: ", round(avg_accuracy * 100, 1), "%"),
            x = "Feature",
            y = "Importance"
          )

        importance_results$plots[[modality]] <- p
      }
    }
  } else {
    # Run on all cells combined
    if (verbose) {
      message("Running feature importance analysis on all cells")
    }

    # Extract cell types
    cell_types <- metadata[, cell_type_var]

    # Run multiple iterations
    iteration_results <- list()
    for (i in 1:n_iterations) {
      iter_result <- runImportanceAnalysis(expr_data, cell_types, method)
      iteration_results[[i]] <- iter_result
    }

    # Combine results across iterations
    combined_importance <- Reduce(`+`, lapply(iteration_results, function(x) x$feature_importance))
    combined_importance <- combined_importance / n_iterations

    # Get top features
    top_features_idx <- order(combined_importance, decreasing = TRUE)[1:min(n_features, length(combined_importance))]
    top_features <- names(combined_importance)[top_features_idx]

    # Calculate stability
    stability_scores <- numeric(length(top_features))
    names(stability_scores) <- top_features

    for (feature in top_features) {
      # Count how many times feature appears in top features across iterations
      feature_ranks <- sapply(iteration_results, function(x) {
        imp <- x$feature_importance
        which(names(imp)[order(imp, decreasing = TRUE)] == feature)
      })

      # Calculate stability as inverse of rank variance
      stability_scores[feature] <- 1 / (1 + stats::var(feature_ranks))
    }

    # Calculate average performance
    avg_accuracy <- mean(sapply(iteration_results, function(x) x$performance$accuracy))

    # Store results
    importance_results$feature_importance[["all"]] <- combined_importance
    importance_results$top_features[["all"]] <- top_features
    importance_results$stability[["all"]] <- stability_scores
    importance_results$performance[["all"]] <- list(accuracy = avg_accuracy)

    # Create visualization
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      # Create data frame for plotting
      plot_data <- data.frame(
        Feature = names(combined_importance[top_features_idx]),
        Importance = combined_importance[top_features_idx],
        Stability = stability_scores[names(combined_importance[top_features_idx])],
        stringsAsFactors = FALSE
      )

      # Create importance plot
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = reorder(Feature, Importance), y = Importance, fill = Stability)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(
          title = "Overall Feature Importance",
          subtitle = paste0("Prediction Accuracy: ", round(avg_accuracy * 100, 1), "%"),
          x = "Feature",
          y = "Importance"
        )

      importance_results$plots[["all"]] <- p
    }
  }

  # Create combined visualization if multiple modalities
  if (modality_specific && length(modalities) > 1 && requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("reshape2", quietly = TRUE)) {
    # Compile top features across modalities
    all_top_features <- unique(unlist(importance_results$top_features))

    # Create feature-modality importance matrix
    feature_modality_matrix <- matrix(0, nrow = length(all_top_features), ncol = length(modalities))
    rownames(feature_modality_matrix) <- all_top_features
    colnames(feature_modality_matrix) <- modalities

    for (modality in modalities) {
      modality_importance <- importance_results$feature_importance[[modality]]
      common_features <- intersect(names(modality_importance), rownames(feature_modality_matrix))
      feature_modality_matrix[common_features, modality] <- modality_importance[common_features]
    }

    # Reshape for plotting
    plot_data <- reshape2::melt(feature_modality_matrix, varnames = c("Feature", "Modality"), value.name = "Importance")

    # Filter to top features
    top_n_overall <- min(20, length(all_top_features))
    top_features_overall <- names(sort(rowSums(feature_modality_matrix), decreasing = TRUE))[1:top_n_overall]
    plot_data_filtered <- plot_data[plot_data$Feature %in% top_features_overall, ]

    # Create heatmap
    heat_plot <- ggplot2::ggplot(plot_data_filtered, ggplot2::aes(x = Modality, y = Feature, fill = Importance)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(
        title = "Feature Importance Across Modalities",
        x = "Modality",
        y = "Feature",
        fill = "Importance"
      )

    importance_results$plots[["comparison"]] <- heat_plot

    # Create overlap plot (UpSet plot)
    if (requireNamespace("UpSetR", quietly = TRUE)) {
      # Create binary presence matrix
      presence_matrix <- matrix(0, nrow = length(all_top_features), ncol = length(modalities))
      rownames(presence_matrix) <- all_top_features
      colnames(presence_matrix) <- modalities

      for (modality in modalities) {
        top_features_mod <- importance_results$top_features[[modality]]
        presence_matrix[top_features_mod, modality] <- 1
      }

      # Convert to data frame for UpSetR
      presence_df <- as.data.frame(presence_matrix)

      # Create UpSet plot
      upset_plot <- UpSetR::upset(presence_df,
                                  nsets = length(modalities),
                                  order.by = "freq",
                                  mainbar.y.label = "Feature Intersection Size",
                                  sets.x.label = "Features Per Modality")

      importance_results$plots[["upset"]] <- upset_plot
    }
  }
  return(importance_results)
}

#' Export or Import Integration Models
#'
#' Provides functions to save and load integration models, allowing users to apply
#' the same integration approach to new datasets.
#'
#' @param integration_results Output from integrateMultiModalData function
#' @param file_path Path where the model will be saved or loaded from
#' @param save_raw_data Logical, whether to save raw data with the model (default: FALSE)
#' @param action Either "save" or "load" (default: "save")
#' @param compression Compression method to use: "gzip", "bzip2", or "xz" (default: "gzip")
#'
#' @return For "save", invisibly returns the file path; for "load", returns the loaded model
#' @export
exportImportIntegrationModel <- function(integration_results = NULL,
                                         file_path,
                                         save_raw_data = FALSE,
                                         action = "save",
                                         compression = "gzip") {

  # Validate parameters
  if (!action %in% c("save", "load")) {
    stop("action must be either 'save' or 'load'")
  }

  if (action == "save") {
    # Check if integration results are provided
    if (is.null(integration_results)) {
      stop("integration_results must be provided for saving")
    }

    # Create export object
    export_obj <- list(
      metadata = integration_results$metadata,
      integration_embeddings = integration_results$integration_embeddings
    )

    # Process integration results based on the object type
    if (inherits(integration_results$integration_results, "Seurat")) {
      # For Seurat objects, extract necessary components
      seurat_obj <- integration_results$integration_results

      # Create a minimal Seurat object
      export_obj$object_type <- "Seurat"

      # Add dimension reductions
      export_obj$reductions <- list()
      for (reduction_name in names(seurat_obj@reductions)) {
        export_obj$reductions[[reduction_name]] <- list(
          embeddings = Seurat::Embeddings(seurat_obj, reduction = reduction_name),
          key = seurat_obj@reductions[[reduction_name]]@key
        )
      }

      # Add metadata
      export_obj$metadata_df <- seurat_obj@meta.data

      # Optionally add raw data
      if (save_raw_data) {
        export_obj$raw_data <- list()
        for (assay_name in names(seurat_obj@assays)) {
          export_obj$raw_data[[assay_name]] <- list(
            counts = as.matrix(Seurat::GetAssayData(seurat_obj, slot = "counts", assay = assay_name)),
            data = as.matrix(Seurat::GetAssayData(seurat_obj, slot = "data", assay = assay_name)),
            scale.data = as.matrix(Seurat::GetAssayData(seurat_obj, slot = "scale.data", assay = assay_name))
          )
        }
      }

    } else if (inherits(integration_results$integration_results, "SingleCellExperiment")) {
      # For SingleCellExperiment objects
      sce_obj <- integration_results$integration_results

      export_obj$object_type <- "SCE"

      # Add reduced dimensions
      export_obj$reductions <- list()
      for (reduction_name in SingleCellExperiment::reducedDimNames(sce_obj)) {
        export_obj$reductions[[reduction_name]] <- SingleCellExperiment::reducedDim(sce_obj, reduction_name)
      }

      # Add metadata
      export_obj$metadata_df <- as.data.frame(SingleCellExperiment::colData(sce_obj))

      # Optionally add raw data
      if (save_raw_data) {
        export_obj$raw_data <- list()
        for (assay_name in SummarizedExperiment::assayNames(sce_obj)) {
          export_obj$raw_data[[assay_name]] <- SummarizedExperiment::assay(sce_obj, assay_name)
        }
      }
    } else {
      warning("Unknown object type. Only saving embeddings and metadata.")
      export_obj$object_type <- "unknown"
    }

    # Save the model
    saveRDS(export_obj, file = file_path, compress = compression)

    # Return file path invisibly
    return(invisible(file_path))

  } else {
    # Load the model
    if (!file.exists(file_path)) {
      stop(paste0("File not found: ", file_path))
    }

    integration_model <- readRDS(file_path)

    # Basic validation
    required_components <- c("metadata", "integration_embeddings", "object_type")
    missing_components <- setdiff(required_components, names(integration_model))

    if (length(missing_components) > 0) {
      warning(paste0("Missing components in the model: ", paste(missing_components, collapse = ", ")))
    }

    return(integration_model)
  }
}

#' Apply Integration Model to New Data
#'
#' Projects new data onto an existing integration model, allowing consistent
#' analysis of new samples within the established integration framework.
#'
#' @param integration_model Output from exportImportIntegrationModel with action="load"
#' @param new_data New data object (must be compatible with the original model)
#' @param reference_batch Name of the reference batch/modality to use for projection
#' @param feature_subset Vector of features to use (default: use all common features)
#' @param scale_new_data Logical, whether to scale new data (default: TRUE)
#' @param method Projection method: "knn" or "linear" (default: "knn")
#' @param n_neighbors Number of neighbors for kNN projection (default: 20)
#' @param seed Random seed for reproducibility (default: 42)
#'
#' @return A list containing projected data and visualization
#' @export
applyIntegrationModel <- function(integration_model,
                                  new_data,
                                  reference_batch,
                                  feature_subset = NULL,
                                  scale_new_data = TRUE,
                                  method = "knn",
                                  n_neighbors = 20,
                                  seed = 42) {

  # Set seed for reproducibility
  set.seed(seed)

  # Check if integration model has necessary components
  required_components <- c("metadata", "integration_embeddings", "object_type")
  missing_components <- setdiff(required_components, names(integration_model))

  if (length(missing_components) > 0) {
    stop(paste0("Missing components in the model: ", paste(missing_components, collapse = ", ")))
  }

  # Check if reference batch exists in the model
  if (!"metadata_df" %in% names(integration_model)) {
    stop("Metadata not found in the integration model")
  }

  if (!"modality" %in% colnames(integration_model$metadata_df)) {
    stop("Modality information not found in the model metadata")
  }

  if (!reference_batch %in% integration_model$metadata_df$modality) {
    stop(paste0("Reference batch '", reference_batch, "' not found in the model"))
  }

  # Get cells from reference batch
  ref_cells <- rownames(integration_model$metadata_df)[integration_model$metadata_df$modality == reference_batch]

  # Extract embeddings for reference cells
  ref_embeddings <- integration_model$integration_embeddings[ref_cells, ]

  # Process new data based on object type
  if (integration_model$object_type == "Seurat") {
    if (!inherits(new_data, "Seurat")) {
      stop("New data must be a Seurat object when the model was built with Seurat")
    }

    # Extract feature space from new data
    expr_matrix <- as.matrix(Seurat::GetAssayData(new_data, slot = "data"))

  } else if (integration_model$object_type == "SCE") {
    if (!inherits(new_data, "SingleCellExperiment")) {
      stop("New data must be a SingleCellExperiment object when the model was built with SCE")
    }

    # Extract feature space from new data
    if ("logcounts" %in% SummarizedExperiment::assayNames(new_data)) {
      expr_matrix <- SummarizedExperiment::assay(new_data, "logcounts")
    } else {
      expr_matrix <- SummarizedExperiment::assay(new_data, 1)
    }

  } else {
    stop(paste0("Unsupported object type: ", integration_model$object_type))
  }

  # If raw data is available in the model, extract reference features
  if ("raw_data" %in% names(integration_model)) {
    # Get reference expression data
    if (integration_model$object_type == "Seurat") {
      assay_name <- names(integration_model$raw_data)[1]  # Use first assay by default
      ref_expr <- integration_model$raw_data[[assay_name]]$data[, ref_cells]

    } else if (integration_model$object_type == "SCE") {
      assay_name <- names(integration_model$raw_data)[1]  # Use first assay by default
      ref_expr <- integration_model$raw_data[[assay_name]][, ref_cells]
    }

    # Find common features
    common_features <- intersect(rownames(ref_expr), rownames(expr_matrix))

    if (length(common_features) < 10) {
      stop("Not enough common features between reference and new data")
    }

    # Subset to common features
    ref_expr <- ref_expr[common_features, , drop = FALSE]
    expr_matrix <- expr_matrix[common_features, , drop = FALSE]

    # Scale data if requested
    if (scale_new_data) {
      # Scale reference data
      ref_expr_scaled <- t(scale(t(ref_expr)))

      # Scale new data
      expr_matrix_scaled <- t(scale(t(expr_matrix)))

      # Replace NaN with 0
      ref_expr_scaled[is.na(ref_expr_scaled)] <- 0
      expr_matrix_scaled[is.na(expr_matrix_scaled)] <- 0

      # Use scaled matrices
      ref_expr <- ref_expr_scaled
      expr_matrix <- expr_matrix_scaled
    }

    # Project new data to the integration space
    if (method == "knn") {
      # Use k-nearest neighbors for projection
      if (requireNamespace("FNN", quietly = TRUE)) {
        # Find k nearest neighbors in feature space
        knn_result <- FNN::get.knnx(
          data = t(ref_expr),
          query = t(expr_matrix),
          k = min(n_neighbors, ncol(ref_expr))
        )

        # Calculate weights based on distances
        weights <- 1 / (knn_result$nn.dist + 1e-10)  # Add small value to avoid division by zero
        weights <- weights / rowSums(weights)  # Normalize weights

        # Project to integration space using weighted average of neighbors
        projected_embeddings <- matrix(0, nrow = nrow(expr_matrix), ncol = ncol(ref_embeddings))
        rownames(projected_embeddings) <- colnames(expr_matrix)
        colnames(projected_embeddings) <- colnames(ref_embeddings)

        for (i in 1:nrow(weights)) {
          neighbors <- knn_result$nn.index[i, ]
          projected_embeddings[i, ] <- colSums(weights[i, ] * ref_embeddings[ref_cells[neighbors], , drop = FALSE])
        }

      } else {
        stop("Package 'FNN' required for kNN projection")
      }
    } else if (method == "linear") {
      # Use linear projection
      if (requireNamespace("irlba", quietly = TRUE)) {
        # Perform PCA on reference data
        pca_result <- irlba::prcomp_irlba(t(ref_expr), n = min(50, ncol(ref_expr) - 1))

        # Project new data onto reference PCA space
        projected_pca <- predict(pca_result, t(expr_matrix))

        # Fit a linear model from PCA to integration space
        fit_models <- list()
        for (j in 1:ncol(ref_embeddings)) {
          fit_models[[j]] <- stats::lm(ref_embeddings[, j] ~ pca_result$x)
        }

        # Project using the linear model
        projected_embeddings <- matrix(0, nrow = nrow(projected_pca), ncol = ncol(ref_embeddings))
        rownames(projected_embeddings) <- colnames(expr_matrix)
        colnames(projected_embeddings) <- colnames(ref_embeddings)

        for (j in 1:ncol(ref_embeddings)) {
          pca_result$x <- projected_pca
          projt <- data.frame(pca_result)
          projected_embeddings[, j] <- predict(fit_models[[j]], newdata = projt)
        }

      } else {
        stop("Package 'irlba' required for linear projection")
      }
    } else {
      stop(paste0("Unsupported projection method: ", method))
    }

    # Add projected data to the new object
    if (integration_model$object_type == "Seurat") {
      # Create a DimReducObject for integration embeddings
      new_data[["integration"]] <- Seurat::CreateDimReducObject(
        embeddings = projected_embeddings,
        key = "integration_"
      )

      # Project to UMAP if available in the model
      if ("reductions" %in% names(integration_model) && "umap" %in% names(integration_model$reductions)) {
        # Get reference UMAP
        ref_umap <- integration_model$reductions$umap$embeddings[ref_cells, ]

        # Project to UMAP using same method
        if (method == "knn") {
          # Use kNN projection
          umap_projected <- matrix(0, nrow = nrow(projected_embeddings), ncol = ncol(ref_umap))
          rownames(umap_projected) <- rownames(projected_embeddings)
          colnames(umap_projected) <- colnames(ref_umap)

          for (i in 1:nrow(weights)) {
            neighbors <- knn_result$nn.index[i, ]
            umap_projected[i, ] <- colSums(weights[i, ] * ref_umap[neighbors, , drop = FALSE])
          }

        } else {
          # Use linear projection
          fit_models_umap <- list()
          for (j in 1:ncol(ref_umap)) {
            fit_models_umap[[j]] <- stats::lm(ref_umap[, j] ~ ref_embeddings)
          }

          umap_projected <- matrix(0, nrow = nrow(projected_embeddings), ncol = ncol(ref_umap))
          rownames(umap_projected) <- rownames(projected_embeddings)
          colnames(umap_projected) <- colnames(ref_umap)

          for (j in 1:ncol(ref_umap)) {
            umap_projected[, j] <- predict(fit_models_umap[[j]],
                                           newdata = data.frame(ref_embeddings = projected_embeddings))
          }
        }

        # Add UMAP to Seurat object
        new_data[["umap"]] <- Seurat::CreateDimReducObject(
          embeddings = umap_projected,
          key = "UMAP_"
        )
      }

    } else if (integration_model$object_type == "SCE") {
      # Add integration embeddings to SCE object
      SingleCellExperiment::reducedDim(new_data, "integration") <- projected_embeddings

      # Project to UMAP if available in the model
      if ("reductions" %in% names(integration_model) && "UMAP" %in% names(integration_model$reductions)) {
        # Get reference UMAP
        ref_umap <- integration_model$reductions$UMAP[ref_cells, ]

        # Project to UMAP using same method
        if (method == "knn") {
          # Use kNN projection
          umap_projected <- matrix(0, nrow = nrow(projected_embeddings), ncol = ncol(ref_umap))
          rownames(umap_projected) <- rownames(projected_embeddings)
          colnames(umap_projected) <- colnames(ref_umap)

          for (i in 1:nrow(weights)) {
            neighbors <- knn_result$nn.index[i, ]
            umap_projected[i, ] <- colSums(weights[i, ] * ref_umap[neighbors, , drop = FALSE])
          }

        } else {
          # Use linear projection
          fit_models_umap <- list()
          for (j in 1:ncol(ref_umap)) {
            fit_models_umap[[j]] <- stats::lm(ref_umap[, j] ~ ref_embeddings)
          }

          umap_projected <- matrix(0, nrow = nrow(projected_embeddings), ncol = ncol(ref_umap))
          rownames(umap_projected) <- rownames(projected_embeddings)
          colnames(umap_projected) <- colnames(ref_umap)

          for (j in 1:ncol(ref_umap)) {
            umap_projected[, j] <- predict(fit_models_umap[[j]],
                                           newdata = data.frame(ref_embeddings = projected_embeddings))
          }
        }

        # Add UMAP to SCE object
        SingleCellExperiment::reducedDim(new_data, "UMAP") <- umap_projected
      }
    }

    # Create results object
    projection_results <- list(
      metadata = list(
        reference_batch = reference_batch,
        method = method,
        n_neighbors = n_neighbors,
        scale_new_data = scale_new_data
      ),
      projected_data = new_data,
      projected_embeddings = projected_embeddings
    )

    # Create visualization
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      # Create data frame for plotting
      if (integration_model$object_type == "Seurat" && "umap" %in% names(new_data@reductions)) {
        # Get UMAP embeddings
        ref_umap <- integration_model$reductions$umap$embeddings
        new_umap <- Seurat::Embeddings(new_data, reduction = "umap")

        plot_df <- data.frame(
          UMAP1 = c(ref_umap[, 1], new_umap[, 1]),
          UMAP2 = c(ref_umap[, 2], new_umap[, 2]),
          Type = c(rep("Reference", nrow(ref_umap)), rep("New", nrow(new_umap))),
          Batch = c(integration_model$metadata_df$modality[match(rownames(ref_umap), rownames(integration_model$metadata_df))],
                    rep("New Data", nrow(new_umap)))
        )

        # Create plot
        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = UMAP1, y = UMAP2, color = Batch, shape = Type)) +
          ggplot2::geom_point(alpha = 0.7, size = 1) +
          ggplot2::theme_minimal() +
          ggplot2::labs(
            title = "Projection of New Data onto Integration Model",
            x = "UMAP1",
            y = "UMAP2"
          )

        projection_results$plot <- p

      } else if (integration_model$object_type == "SCE" && "UMAP" %in% SingleCellExperiment::reducedDimNames(new_data)) {
        # Get UMAP embeddings
        ref_umap <- integration_model$reductions$UMAP
        new_umap <- SingleCellExperiment::reducedDim(new_data, "UMAP")

        plot_df <- data.frame(
          UMAP1 = c(ref_umap[, 1], new_umap[, 1]),
          UMAP2 = c(ref_umap[, 2], new_umap[, 2]),
          Type = c(rep("Reference", nrow(ref_umap)), rep("New", nrow(new_umap))),
          Batch = c(integration_model$metadata_df$modality[match(rownames(ref_umap), rownames(integration_model$metadata_df))],
                    rep("New Data", nrow(new_umap)))
        )

        # Create plot
        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = UMAP1, y = UMAP2, color = Batch, shape = Type)) +
          ggplot2::geom_point(alpha = 0.7, size = 1) +
          ggplot2::theme_minimal() +
          ggplot2::labs(
            title = "Projection of New Data onto Integration Model",
            x = "UMAP1",
            y = "UMAP2"
          )

        projection_results$plot <- p
      }
    }

    return(projection_results)

  } else {
    stop("Raw data not available in the model. Cannot project new data.")
  }
}
