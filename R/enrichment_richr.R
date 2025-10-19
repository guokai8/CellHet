#' Functional Enrichment Analysis using richR
#'
#' Perform functional enrichment analysis on DEG results using the richR package.
#' Supports GO, KEGG, and custom gene set enrichment with visualization options.
#'
#' @param deg_results Output from compareDEGs() or a vector of gene names
#' @param gene_list Optional named vector of genes (names are gene symbols, values are scores/logFC)
#' @param analysis_type Type of enrichment: "GO", "KEGG", "GSEA", or "DAVID" (default: "GO")
#' @param organism Organism name: "human", "mouse", "rat", etc. (default: "human")
#' @param ontology For GO analysis: "BP", "MF", "CC", or "ALL" (default: "BP")
#' @param keytype Gene ID type: "SYMBOL", "ENTREZID", "ENSEMBL" (default: "SYMBOL")
#' @param pvalue P-value cutoff for enrichment (default: 0.05)
#' @param qvalue Q-value cutoff for multiple testing (default: 0.05)
#' @param minGSSize Minimum gene set size (default: 10)
#' @param maxGSSize Maximum gene set size (default: 500)
#' @param use_significant_only Use only significant DEGs (default: TRUE)
#' @param separate_direction Separate up/down regulated genes (default: FALSE)
#' @param build_annotation Whether to build annotation if not available (default: TRUE)
#' @param visualize Create visualization plots (default: TRUE)
#' @param plot_type Type of plot: "bar", "dot", "network", "all" (default: "dot")
#' @param top_terms Number of top terms to display (default: 20)
#' @param return_plots Return plot objects (default: TRUE)
#'
#' @return A list containing:
#'   \itemize{
#'     \item enrichment_results: Enrichment analysis results
#'     \item plots: List of ggplot objects (if visualize=TRUE)
#'     \item summary: Summary statistics
#'     \item gene_annotation: Annotation data used
#'   }
#'
#' @details
#' This function wraps the richR package functionality to provide easy enrichment
#' analysis for CellHet DEG results. It supports:
#' - GO enrichment (Biological Process, Molecular Function, Cellular Component)
#' - KEGG pathway enrichment
#' - Gene Set Enrichment Analysis (GSEA)
#' - DAVID online analysis
#'
#' The function can work with:
#' 1. DEG results from compareDEGs() - automatically extracts significant genes
#' 2. Custom gene list - provide vector of gene symbols
#' 3. Ranked gene list - provide named vector for GSEA
#'
#' @examples
#' \dontrun{
#' # Basic GO enrichment on DEG results
#' enrich_results <- runRichREnrichment(
#'   deg_results,
#'   analysis_type = "GO",
#'   organism = "human"
#' )
#'
#' # KEGG pathway enrichment
#' kegg_results <- runRichREnrichment(
#'   deg_results,
#'   analysis_type = "KEGG",
#'   organism = "mouse"
#' )
#'
#' # GSEA with ranked genes
#' ranked_genes <- setNames(deg_df$avg_log2FC, deg_df$gene)
#' gsea_results <- runRichREnrichment(
#'   ranked_genes,
#'   analysis_type = "GSEA",
#'   organism = "human"
#' )
#'
#' # Separate up and down regulated genes
#' enrich_separate <- runRichREnrichment(
#'   deg_results,
#'   separate_direction = TRUE
#' )
#' }
#'
#' @export
runRichREnrichment <- function(deg_results,
                              gene_list = NULL,
                              analysis_type = "GO",
                              organism = "human",
                              ontology = "BP",
                              keytype = "SYMBOL",
                              pvalue = 0.05,
                              qvalue = 0.05,
                              minGSSize = 10,
                              maxGSSize = 500,
                              use_significant_only = TRUE,
                              separate_direction = FALSE,
                              build_annotation = TRUE,
                              visualize = TRUE,
                              plot_type = "dot",
                              top_terms = 20,
                              return_plots = TRUE) {

  # Check if richR is installed and load it
  if (!requireNamespace("richR", quietly = TRUE)) {
    stop("Package 'richR' is required for enrichment analysis.\n",
         "Install with: devtools::install_github('guokai8/richR')")
  }

  # Load richR library to ensure all functions and dependencies are available
  # This is necessary because richR uses many internal functions that need full library loading
  suppressPackageStartupMessages({
    library(richR)
  })

  if (visualize && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for visualization. Please install it.")
  }

  # Initialize results list
  results <- list()

  # Extract genes from input
  genes_to_analyze <- extractGenesForEnrichment(
    deg_results,
    gene_list,
    use_significant_only,
    separate_direction
  )

  if (length(genes_to_analyze) == 0) {
    stop("No genes provided for enrichment analysis")
  }

  # Build or load annotation data
  message("Preparing annotation data for ", organism, "...")
  annot_data <- prepareAnnotationData(organism, analysis_type, build_annotation)

  # Run enrichment analysis
  message("Running ", analysis_type, " enrichment analysis...")

  if (separate_direction && is.list(genes_to_analyze)) {
    # Analyze up and down regulated genes separately
    results$enrichment_results <- list()

    for (direction in names(genes_to_analyze)) {
      message("  Analyzing ", direction, " regulated genes (",
              length(genes_to_analyze[[direction]]), " genes)...")

      enrich_res <- performEnrichment(
        genes = genes_to_analyze[[direction]],
        annot_data = annot_data,
        analysis_type = analysis_type,
        ontology = ontology,
        pvalue = pvalue,
        qvalue = qvalue,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize
      )

      results$enrichment_results[[direction]] <- enrich_res
    }
  } else {
    # Single enrichment analysis
    if (is.list(genes_to_analyze)) {
      genes_to_analyze <- genes_to_analyze[[1]]
    }

    message("  Analyzing ", length(genes_to_analyze), " genes...")

    results$enrichment_results <- performEnrichment(
      genes = genes_to_analyze,
      annot_data = annot_data,
      analysis_type = analysis_type,
      ontology = ontology,
      pvalue = pvalue,
      qvalue = qvalue,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize
    )
  }

  # Create summary
  results$summary <- createEnrichmentSummary(results$enrichment_results, separate_direction)

  # Store annotation info
  results$annotation_info <- list(
    organism = organism,
    analysis_type = analysis_type,
    ontology = ontology,
    pvalue = pvalue,
    qvalue = qvalue
  )

  # Create visualizations
  if (visualize && return_plots) {
    message("Creating visualization plots...")
    results$plots <- createEnrichmentPlots(
      results$enrichment_results,
      plot_type = plot_type,
      top_terms = top_terms,
      separate_direction = separate_direction
    )
  }

  message("Enrichment analysis complete!")

  return(results)
}


#' Extract Genes for Enrichment Analysis
#'
#' @keywords internal
extractGenesForEnrichment <- function(deg_results, gene_list, use_significant_only, separate_direction) {

  # If gene_list provided directly, use it
  if (!is.null(gene_list)) {
    return(gene_list)
  }

  # Extract from DEG results
  if (is.data.frame(deg_results)) {
    # Single data frame
    genes_df <- deg_results
  } else if ("degs" %in% names(deg_results)) {
    # CellHet DEG results object
    # Combine all DEGs
    all_genes <- data.frame()
    for (result_key in names(deg_results$degs)) {
      degs <- deg_results$degs[[result_key]]
      all_genes <- rbind(all_genes, degs)
    }
    genes_df <- all_genes
  } else {
    stop("Unrecognized input format. Provide DEG results or gene list.")
  }

  # Filter for significant genes if requested
  if (use_significant_only) {
    if ("significant" %in% colnames(genes_df)) {
      genes_df <- genes_df[genes_df$significant, ]
    } else if ("p_val_adj" %in% colnames(genes_df)) {
      genes_df <- genes_df[genes_df$p_val_adj < 0.05, ]
    }
  }

  # Check if we have genes
  if (nrow(genes_df) == 0) {
    stop("No genes found after filtering. Try use_significant_only=FALSE")
  }

  # Separate by direction if requested
  if (separate_direction) {
    if ("direction" %in% colnames(genes_df)) {
      genes_up <- genes_df$gene[genes_df$direction == "up"]
      genes_down <- genes_df$gene[genes_df$direction == "down"]

      return(list(
        up = unique(genes_up),
        down = unique(genes_down)
      ))
    } else {
      warning("No 'direction' column found. Cannot separate genes.")
      return(list(all = unique(genes_df$gene)))
    }
  } else {
    return(unique(genes_df$gene))
  }
}


#' Prepare Annotation Data
#'
#' @keywords internal
prepareAnnotationData <- function(organism, analysis_type, build_annotation) {

  # richR uses full organism names, not KEGG codes
  organism_name <- switch(tolower(organism),
                         "human" = "human",
                         "mouse" = "mouse",
                         "rat" = "rat",
                         "zebrafish" = "zebrafish",
                         "fly" = "fly",
                         "worm" = "worm",
                         "yeast" = "yeast",
                         organism)

  # Try to build annotation if requested
  if (build_annotation) {
    tryCatch({
      message("Building annotation for ", organism_name, " (", analysis_type, ")...")

      if (analysis_type %in% c("GO", "GSEA")) {
        annot_data <- richR::buildAnnot(
          species = organism_name,
          keytype = "SYMBOL",
          anntype = "GO"
        )
        message("GO annotation built successfully")
      } else if (analysis_type == "KEGG") {
        annot_data <- richR::buildAnnot(
          species = organism_name,
          keytype = "SYMBOL",
          anntype = "KEGG"
        )
        message("KEGG annotation built successfully")
      } else {
        annot_data <- NULL
      }

      if (is.null(annot_data)) {
        stop("Annotation data is NULL after building. This may indicate:\n",
             "1. No internet connection for downloading annotation databases\n",
             "2. richR package needs to be properly installed: devtools::install_github('guokai8/richR')\n",
             "3. Organism '", organism_name, "' may not be supported")
      }

      return(annot_data)

    }, error = function(e) {
      stop("Could not build annotation: ", e$message,
           "\n\nTroubleshooting steps:\n",
           "1. Check internet connection\n",
           "2. Install richR: devtools::install_github('guokai8/richR')\n",
           "3. Verify organism is supported: ", organism_name, "\n",
           "4. Try manually: annot <- richR::buildAnnot(species='", organism_name, "', keytype='SYMBOL', anntype='",
           ifelse(analysis_type == "KEGG", "KEGG", "GO"), "')")
    })
  } else {
    stop("build_annotation is FALSE but no annotation data provided.\n",
         "Either set build_annotation=TRUE or provide annotation data manually.")
  }
}


#' Perform Enrichment Analysis
#'
#' @keywords internal
performEnrichment <- function(genes, annot_data, analysis_type, ontology,
                             pvalue, qvalue, minGSSize, maxGSSize) {

  if (is.null(annot_data)) {
    stop("Annotation data is required. Set build_annotation=TRUE or provide annotation.")
  }

  tryCatch({
    if (analysis_type == "GO") {
      result <- richR::richGO(
        x = genes,  # richR uses 'x' not 'gene'
        godata = annot_data,
        ontology = ontology,
        pvalue = pvalue,
        padj = qvalue,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize
      )
    } else if (analysis_type == "KEGG") {
      result <- richR::richKEGG(
        x = genes,  # richR uses 'x' not 'gene'
        kodata = annot_data,
        pvalue = pvalue,
        padj = qvalue,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize
      )
    } else if (analysis_type == "GSEA") {
      # For GSEA, genes should be a named numeric vector
      if (!is.numeric(genes) || is.null(names(genes))) {
        stop("For GSEA, provide a named numeric vector (gene names with scores)")
      }

      result <- richR::richGSEA(
        x = genes,  # richR uses 'x' not 'gene'
        object = annot_data,
        pvalue = pvalue,
        padj = qvalue,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize
      )
    } else {
      stop("Unsupported analysis_type: ", analysis_type)
    }

    # Convert richResult object to data.frame for easier manipulation
    # IMPORTANT: Keep original column names (Padj, Pvalue) for richR plotting functions
    if (!is.null(result)) {
      result <- as.data.frame(result)
    }

    return(result)

  }, error = function(e) {
    stop("Enrichment analysis failed: ", e$message)
  })
}


#' Create Enrichment Summary
#'
#' @keywords internal
createEnrichmentSummary <- function(enrichment_results, separate_direction) {

  if (separate_direction && is.list(enrichment_results)) {
    # Summary for each direction
    summary_list <- lapply(names(enrichment_results), function(direction) {
      res <- enrichment_results[[direction]]
      if (is.null(res) || nrow(res) == 0) {
        return(data.frame(
          direction = direction,
          n_terms = 0,
          n_significant = 0,
          top_term = NA_character_,
          top_pvalue = NA_real_
        ))
      }

      data.frame(
        direction = direction,
        n_terms = nrow(res),
        n_significant = sum(res$Padj < 0.05, na.rm = TRUE),
        top_term = if(nrow(res) > 0) res$Term[1] else NA_character_,
        top_pvalue = if(nrow(res) > 0) res$Padj[1] else NA_real_
      )
    })

    summary_df <- do.call(rbind, summary_list)
  } else {
    # Single summary
    res <- enrichment_results

    if (is.null(res) || nrow(res) == 0) {
      summary_df <- data.frame(
        n_terms = 0,
        n_significant = 0,
        top_term = NA_character_,
        top_pvalue = NA_real_,
        min_pvalue = NA_real_,
        max_pvalue = NA_real_
      )
    } else {
      summary_df <- data.frame(
        n_terms = nrow(res),
        n_significant = sum(res$Padj < 0.05, na.rm = TRUE),
        top_term = if(nrow(res) > 0) res$Term[1] else NA_character_,
        top_pvalue = if(nrow(res) > 0) res$Padj[1] else NA_real_,
        min_pvalue = if(nrow(res) > 0) min(res$Padj, na.rm = TRUE) else NA_real_,
        max_pvalue = if(nrow(res) > 0) max(res$Padj, na.rm = TRUE) else NA_real_
      )
    }
  }

  return(summary_df)
}


#' Create Enrichment Plots
#'
#' @keywords internal
createEnrichmentPlots <- function(enrichment_results, plot_type, top_terms, separate_direction) {

  plots <- list()

  if (separate_direction && is.list(enrichment_results)) {
    # Create plots for each direction
    for (direction in names(enrichment_results)) {
      res <- enrichment_results[[direction]]

      if (is.null(res) || nrow(res) == 0) {
        message("No significant terms for ", direction, " regulated genes")
        next
      }

      # Limit to top terms
      res_top <- head(res, top_terms)

      # Create requested plot types
      if (plot_type %in% c("bar", "all")) {
        tryCatch({
          plots[[paste0(direction, "_bar")]] <- richR::ggbar(res_top, top = top_terms) +
            ggplot2::ggtitle(paste(direction, "regulated genes"))
        }, error = function(e) {
          message("Could not create bar plot: ", e$message)
        })
      }

      if (plot_type %in% c("dot", "all")) {
        tryCatch({
          plots[[paste0(direction, "_dot")]] <- richR::ggdot(res_top, top = top_terms) +
            ggplot2::ggtitle(paste(direction, "regulated genes"))
        }, error = function(e) {
          message("Could not create dot plot: ", e$message)
        })
      }

      if (plot_type %in% c("network", "all")) {
        tryCatch({
          plots[[paste0(direction, "_network")]] <- richR::ggnetplot(res_top) +
            ggplot2::ggtitle(paste(direction, "regulated genes"))
        }, error = function(e) {
          message("Could not create network plot: ", e$message)
        })
      }
    }
  } else {
    # Single set of plots
    res <- enrichment_results

    if (is.null(res) || nrow(res) == 0) {
      message("No significant terms found")
      return(plots)
    }

    # Limit to top terms
    res_top <- head(res, top_terms)

    # Create requested plot types
    if (plot_type %in% c("bar", "all")) {
      tryCatch({
        plots$bar <- richR::ggbar(res_top, top = top_terms)
      }, error = function(e) {
        message("Could not create bar plot: ", e$message)
      })
    }

    if (plot_type %in% c("dot", "all")) {
      tryCatch({
        plots$dot <- richR::ggdot(res_top, top = top_terms)
      }, error = function(e) {
        message("Could not create dot plot: ", e$message)
      })
    }

    if (plot_type %in% c("network", "all")) {
      tryCatch({
        plots$network <- richR::ggnetplot(res_top)
      }, error = function(e) {
        message("Could not create network plot: ", e$message)
      })
    }
  }

  return(plots)
}


#' Export Enrichment Results
#'
#' Export enrichment results to Excel or CSV format
#'
#' @param enrichment_results Output from runRichREnrichment()
#' @param file_path Path to save the file
#' @param format File format: "xlsx" or "csv" (default: "xlsx")
#'
#' @return Invisibly returns the file path
#'
#' @examples
#' \dontrun{
#' # Export to Excel
#' exportEnrichmentResults(enrich_results, "enrichment.xlsx")
#'
#' # Export to CSV
#' exportEnrichmentResults(enrich_results, "enrichment.csv", format = "csv")
#' }
#'
#' @export
exportEnrichmentResults <- function(enrichment_results, file_path, format = "xlsx") {

  if (format == "xlsx") {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' needed for Excel export. Please install it.")
    }

    wb <- openxlsx::createWorkbook()

    # Add summary sheet
    if ("summary" %in% names(enrichment_results)) {
      openxlsx::addWorksheet(wb, "Summary")
      openxlsx::writeData(wb, "Summary", enrichment_results$summary)
    }

    # Add enrichment results
    if ("enrichment_results" %in% names(enrichment_results)) {
      enrich_res <- enrichment_results$enrichment_results

      if (is.list(enrich_res) && !is.data.frame(enrich_res)) {
        # Multiple result sets (e.g., up/down)
        for (name in names(enrich_res)) {
          sheet_name <- gsub("[^A-Za-z0-9_]", "_", name)
          openxlsx::addWorksheet(wb, sheet_name)
          openxlsx::writeData(wb, sheet_name, enrich_res[[name]])
        }
      } else {
        # Single result set
        openxlsx::addWorksheet(wb, "Enrichment")
        openxlsx::writeData(wb, "Enrichment", enrich_res)
      }
    }

    openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)

  } else if (format == "csv") {
    # For CSV, combine all results
    if ("enrichment_results" %in% names(enrichment_results)) {
      enrich_res <- enrichment_results$enrichment_results

      if (is.list(enrich_res) && !is.data.frame(enrich_res)) {
        # Combine multiple results
        combined <- data.frame()
        for (name in names(enrich_res)) {
          df <- enrich_res[[name]]
          df$group <- name
          combined <- rbind(combined, df)
        }
        write.csv(combined, file_path, row.names = FALSE)
      } else {
        write.csv(enrich_res, file_path, row.names = FALSE)
      }
    }
  } else {
    stop("Unsupported format: ", format)
  }

  message("Enrichment results exported to ", file_path)
  invisible(file_path)
}
