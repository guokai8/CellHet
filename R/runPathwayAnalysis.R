#' Perform Pathway Enrichment Analysis on DEGs using richR
#'
#' @param deg_results Output from compareDEGs function or findDifferentialGenes function
#' @param annot_data Annotation data created by buildAnnot from richR
#' @param annot_type Type of annotation: "GO", "KEGG", "KEGGM", or "MSIGDB"
#' @param ontology For GO analysis, the ontology to use: "BP", "MF", or "CC" (default: "BP")
#' @param pvalue P-value threshold for pathway significance (default: 0.05)
#' @param padj Adjusted p-value threshold for pathway significance (default: 0.05)
#' @param padj.method Method for p-value adjustment (default: "BH")
#' @param direction Filter for gene direction: "up", "down", or "both" (default: "both")
#' @param min_genes Minimum number of genes required in a pathway (default: 5)
#' @param max_genes Maximum number of genes allowed in a pathway (default: 500)
#' @param use_gsea Whether to use Gene Set Enrichment Analysis instead of ORA (default: FALSE)
#' @param verbose Whether to show progress messages (default: TRUE)
#'
#' @return A list with pathway analysis results for each comparison
#' @export
runPathwayAnalysis <- function(deg_results,
                               annot_data,
                               annot_type = c("GO", "KEGG", "KEGGM", "MSIGDB"),
                               ontology = "BP",
                               pvalue = 0.05,
                               padj = 0.05,
                               padj.method = "BH",
                               direction = "both",
                               min_genes = 5,
                               max_genes = 500,
                               use_gsea = FALSE,
                               verbose = TRUE) {

  # Check if richR is installed
  if (!requireNamespace("richR", quietly = TRUE)) {
    stop("Package 'richR' is required for pathway analysis. Please install it with:\n",
         "devtools::install_github('guokai8/richR')")
  }

  # Match annotation type
  annot_type <- match.arg(annot_type)

  # Extract data based on input type
  if (inherits(deg_results, "list")) {
    if ("all_degs" %in% names(deg_results)) {
      # Data from findDifferentialGenes
      all_degs <- deg_results$all_degs
    } else if ("degs" %in% names(deg_results)) {
      # Data from compareDEGs
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

        # Split comparison into groups if needed
        if (!("group1" %in% colnames(degs)) && !("group2" %in% colnames(degs))) {
          comp_parts <- strsplit(comparison, "_vs_")[[1]]
          if (length(comp_parts) == 2) {
            degs$group1 <- comp_parts[1]
            degs$group2 <- comp_parts[2]
          }
        }

        # Append to combined data frame
        all_degs <- rbind(all_degs, degs)
      }
    } else if (length(deg_results) > 0 && all(sapply(deg_results, function(x) is.data.frame(x)))) {
      # List of data frames - compile them
      all_degs <- do.call(rbind, deg_results)
    } else {
      stop("Unrecognized input format. Must be output from compareDEGs/findDifferentialGenes or a list of data frames.")
    }
  } else if (inherits(deg_results, "data.frame")) {
    # Input is already a data frame
    all_degs <- deg_results
  } else {
    stop("Input must be a data frame or a list containing DEG results")
  }

  # Check required columns
  req_cols <- c("gene", "significant")
  if (!all(req_cols %in% colnames(all_degs))) {
    stop("DEG results must contain the columns: ", paste(req_cols, collapse = ", "))
  }

  # For GSEA, we need avg_log2FC values for all genes
  if (use_gsea && (!("avg_log2FC" %in% colnames(all_degs)))) {
    stop("For GSEA, DEG results must contain the column: avg_log2FC")
  }

  # Check if there's a proper direction column
  if (direction != "both") {
    if ("direction" %in% colnames(all_degs)) {
      # Use direction column
    } else if ("avg_log2FC" %in% colnames(all_degs)) {
      # Determine direction from log fold change
      all_degs$direction <- ifelse(all_degs$avg_log2FC > 0, "up", "down")
    } else {
      warning("Cannot filter by direction - no direction or avg_log2FC column found. Using all significant genes.")
      direction <- "both"
    }
  }

  # Check for passes_lfc_threshold column
  has_lfc_threshold_info <- "passes_lfc_threshold" %in% colnames(all_degs)

  # Get unique cell type and comparison combinations
  if ("cell_type" %in% colnames(all_degs) && "comparison" %in% colnames(all_degs)) {
    ct_comp_pairs <- unique(all_degs[, c("cell_type", "comparison")])
  } else if ("cell_type" %in% colnames(all_degs)) {
    ct_comp_pairs <- data.frame(
      cell_type = unique(all_degs$cell_type),
      comparison = "all"
    )
  } else if ("comparison" %in% colnames(all_degs)) {
    ct_comp_pairs <- data.frame(
      cell_type = "all",
      comparison = unique(all_degs$comparison)
    )
  } else {
    ct_comp_pairs <- data.frame(
      cell_type = "all",
      comparison = "all"
    )
  }

  # Initialize results storage
  pathway_results <- list()

  # Counter for progress reporting
  total_combos <- nrow(ct_comp_pairs)
  counter <- 0

  # Process each cell type and comparison
  for (i in 1:nrow(ct_comp_pairs)) {
    ct <- ct_comp_pairs$cell_type[i]
    comp <- ct_comp_pairs$comparison[i]

    counter <- counter + 1
    if (verbose) {
      message(sprintf("Processing %s - %s (%d of %d)", ct, comp, counter, total_combos))
    }

    # Filter DEGs for this cell type and comparison
    if (ct == "all" && comp == "all") {
      subset_degs <- all_degs
    } else if (ct == "all") {
      subset_degs <- all_degs[all_degs$comparison == comp, ]
    } else if (comp == "all") {
      subset_degs <- all_degs[all_degs$cell_type == ct, ]
    } else {
      subset_degs <- all_degs[all_degs$cell_type == ct & all_degs$comparison == comp, ]
    }

    if (!use_gsea) {
      # For ORA, filter for significant genes that pass thresholds
      sig_degs <- subset_degs[subset_degs$significant, ]

      # Apply logFC threshold if available
      if (has_lfc_threshold_info) {
        sig_degs <- sig_degs[sig_degs$passes_lfc_threshold, ]
      }

      # Filter by direction if specified
      if (direction == "up") {
        sig_degs <- sig_degs[sig_degs$direction == "up", ]
      } else if (direction == "down") {
        sig_degs <- sig_degs[sig_degs$direction == "down", ]
      }

      # Extract gene list
      genes <- unique(sig_degs$gene)

      # Skip if no significant genes found
      if (length(genes) == 0) {
        if (verbose) {
          message(sprintf("No significant genes found for %s - %s. Skipping.", ct, comp))
        }
        next
      }
    } else {
      # For GSEA, prepare a ranked gene list using all genes, not just significant ones
      # Filter by direction if specified, but keep all genes within that direction
      # For GSEA, ensure we have enough genes
      if (nrow(subset_degs) < 100) {
        warning(sprintf("Only %d genes available for GSEA for %s - %s. Results may be unreliable.",
                        nrow(subset_degs), ct, comp))
        if (nrow(subset_degs) == 0) {
          next  # Skip if no genes
        }
      }

      # Prepare gene ranking for GSEA
      gene_ranks <- subset_degs$avg_log2FC
      names(gene_ranks) <- subset_degs$gene
    }

    # Create result key
    result_key <- paste(ct, comp, sep = "__")

    # Run appropriate enrichment analysis based on annotation type
    tryCatch({
      if (use_gsea) {
        # For GSEA, we need to run richGSEA with the ranked gene list
        if (annot_type == "GO") {
          enrichment_result <- richR::richGSEA(
            x = gene_ranks,
            object = annot_data,
            ontology = ontology,
            pvalue = pvalue,
            padj.method = padj.method,
            minSize = min_genes,
            maxSize = max_genes
          )
        } else if (annot_type == "KEGG" || annot_type == "KEGGM") {
          enrichment_result <- richR::richGSEA(
            x = gene_ranks,
            object = annot_data,
            pvalue = pvalue,
            padj.method = padj.method,
            minSize = min_genes,
            maxSize = max_genes
          )
        } else if (annot_type == "MSIGDB") {
          enrichment_result <- richR::richGSEA(
            x = gene_ranks,
            object = annot_data,
            pvalue = pvalue,
            padj.method = padj.method,
            minSize = min_genes,
            maxSize = max_genes
          )
        }
      } else {
        # Run over-representation analysis
        if (annot_type == "GO") {
          enrichment_result <- richR::richGO(
            x = genes,
            godata = annot_data,
            ontology = ontology,
            pvalue = pvalue,
            padj.method = padj.method,
            minGSSize = min_genes,
            maxGSSize = max_genes
          )
        } else if (annot_type == "KEGG") {
          enrichment_result <- richR::richKEGG(
            x = genes,
            kodata = annot_data,
            pvalue = pvalue,
            padj.method = padj.method,
            minGSSize = min_genes,
            maxGSSize = max_genes
          )
        } else if (annot_type == "KEGGM") {
          # For KEGGM, use enrich function
          enrichment_result <- richR::enrich(
            x = genes,
            object = annot_data,
            pvalue = pvalue,
            padj.method = padj.method,
            minGSSize = min_genes,
            maxGSSize = max_genes
          )
        } else if (annot_type == "MSIGDB") {
          # For MSIGDB, use enrich function
          enrichment_result <- richR::enrich(
            x = genes,
            object = annot_data,
            pvalue = pvalue,
            padj.method = padj.method,
            minGSSize = min_genes,
            maxGSSize = max_genes
          )
        }
      }

      # Apply p-value and padj filters
      if (!is.null(enrichment_result) && nrow(enrichment_result) > 0) {
        # Convert enrichment result to data frame if it's not already
        enrichment_df <- as.data.frame(enrichment_result)

        # Apply p-value thresholds - check for both Pvalue and pvalue column names
        if ("Pvalue" %in% colnames(enrichment_df)) {
          enrichment_df <- enrichment_df[enrichment_df$Pvalue <= pvalue, ]
        } else if ("pvalue" %in% colnames(enrichment_df)) {
          enrichment_df <- enrichment_df[enrichment_df$pvalue <= pvalue, ]
        }

        # Apply adjusted p-value thresholds - check for both Padj and p.adjust column names
        if ("Padj" %in% colnames(enrichment_df)) {
          enrichment_df <- enrichment_df[enrichment_df$Padj <= padj, ]
        } else if ("p.adjust" %in% colnames(enrichment_df)) {
          enrichment_df <- enrichment_df[enrichment_df$p.adjust <= padj, ]
        }

        # Skip if no results pass filters
        if (nrow(enrichment_df) == 0) {
          if (verbose) {
            message(sprintf("No significant pathways found for %s - %s after filtering.", ct, comp))
          }
          next
        }

        # Store results - include ALL genes in the genes list for GSEA
        pathway_results[[result_key]] <- list(
          enrichment = enrichment_df,
          genes = if(use_gsea) names(gene_ranks) else genes,
          num_genes = if(use_gsea) length(gene_ranks) else length(genes),
          direction = direction,
          cell_type = ct,
          comparison = comp,
          annot_type = annot_type,
          method = if(use_gsea) "GSEA" else "ORA"
        )

        # Get gene-pathway details if available
        if (!use_gsea) {
          tryCatch({
            details <- richR::detail(enrichment_result)
            pathway_results[[result_key]]$details <- details
          }, error = function(e) {
            warning("Error getting pathway details for ", result_key, ": ", e$message)
          })
        }
      } else {
        if (verbose) {
          message(sprintf("No significant pathways found for %s - %s.", ct, comp))
        }
      }
    }, error = function(e) {
      warning("Error in pathway analysis for ", result_key, ": ", e$message)
    })
  }

  # Return results
  return(list(
    pathway_results = pathway_results,
    metadata = list(
      annot_type = annot_type,
      ontology = if(annot_type == "GO") ontology else NULL,
      direction = direction,
      method = if(use_gsea) "GSEA" else "ORA",
      pvalue = pvalue,
      padj = padj,
      min_genes = min_genes,
      max_genes = max_genes
    )
  ))
}
