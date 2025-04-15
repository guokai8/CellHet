# runPathwayAnalysis.R

#' Run Simple Gene Set Enrichment Analysis
#'
#' Performs over-representation analysis on gene lists using
#' hypergeometric test without requiring clusterProfiler package.
#'
#' @param gene_list A vector of gene symbols to test for enrichment
#' @param gene_sets A named list where each element is a vector of genes in a gene set
#' @param universe A vector of all possible background genes (default: union of all gene sets)
#' @param min_genes Minimum number of genes required for a pathway (default: 5)
#' @param max_genes Maximum number of genes allowed for a pathway (default: 500)
#' @param p_adj_method Method for p-value adjustment (default: "BH")
#' @param p_val_cutoff P-value threshold for significance (default: 0.05)
#'
#' @return A data frame with enrichment results
#' @keywords internal
simple_ora <- function(gene_list,
                       gene_sets,
                       universe = NULL,
                       min_genes = 5,
                       max_genes = 500,
                       p_adj_method = "BH",
                       p_val_cutoff = 0.05) {

  # Validate input
  if (length(gene_list) == 0) {
    return(data.frame())
  }

  # Create a universe of genes if not provided
  if (is.null(universe)) {
    universe <- unique(unlist(gene_sets))
  }

  # Filter for genes in universe to ensure valid analysis
  gene_list <- gene_list[gene_list %in% universe]

  if (length(gene_list) == 0) {
    warning("No genes in the input list are found in the universe")
    return(data.frame())
  }

  # Initialize results table
  results <- data.frame(
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer(),
    stringsAsFactors = FALSE
  )

  # Size of gene list and universe
  n_gene_list <- length(gene_list)
  n_universe <- length(universe)

  # Loop through each gene set
  for (set_name in names(gene_sets)) {
    set_genes <- gene_sets[[set_name]]

    # Filter for genes in universe
    set_genes <- set_genes[set_genes %in% universe]

    # Skip if gene set doesn't meet size criteria
    if (length(set_genes) < min_genes || length(set_genes) > max_genes) {
      next
    }

    # Get overlapping genes
    overlap_genes <- intersect(gene_list, set_genes)
    n_overlap <- length(overlap_genes)

    # Skip if no overlapping genes
    if (n_overlap == 0) {
      next
    }

    # Number of genes in the gene set
    n_gene_set <- length(set_genes)

    # Hypergeometric test
    # phyper(q, m, n, k, lower.tail = FALSE) computes P(X > q)
    # q = number of overlapping genes - 1 (for > q)
    # m = number of genes in gene set
    # n = number of genes in universe NOT in gene set
    # k = number of genes in gene list
    p_value <- stats::phyper(
      q = n_overlap - 1,
      m = n_gene_set,
      n = n_universe - n_gene_set,
      k = n_gene_list,
      lower.tail = FALSE
    )

    # Add to results table
    results <- rbind(results, data.frame(
      ID = set_name,
      Description = set_name,  # Can be replaced with actual descriptions if available
      GeneRatio = paste0(n_overlap, "/", n_gene_list),
      BgRatio = paste0(n_gene_set, "/", n_universe),
      pvalue = p_value,
      p.adjust = NA,  # Will adjust later
      qvalue = NA,    # Will set equal to p.adjust
      geneID = paste(overlap_genes, collapse = "/"),
      Count = n_overlap,
      stringsAsFactors = FALSE
    ))
  }

  # If no enriched gene sets found, return empty data frame
  if (nrow(results) == 0) {
    return(results)
  }

  # Adjust p-values
  results$p.adjust <- stats::p.adjust(results$pvalue, method = p_adj_method)
  results$qvalue <- results$p.adjust

  # Sort by p-value
  results <- results[order(results$pvalue), ]

  # Filter by adjusted p-value
  results <- results[results$p.adjust < p_val_cutoff, ]

  return(results)
}

#' Perform Pathway Enrichment Analysis on DEGs
#'
#' @param deg_results Output from compareDEGs function
#' @param gene_sets Named list of gene sets for enrichment analysis
#' @param background Optional vector of background genes
#' @param min_genes Minimum number of genes required for a pathway
#' @param p_adj_cutoff Adjusted p-value threshold for pathway significance
#' @param direction Filter for gene direction ("up", "down", "both")
#'
#' @return A list with pathway analysis results for each comparison
#' @export
runPathwayAnalysis <- function(deg_results,
                               gene_sets,
                               background = NULL,
                               min_genes = 5,
                               p_adj_cutoff = 0.05,
                               direction = "both") {

  # Extract DEG results
  deg_list <- deg_results$degs

  if (length(deg_list) == 0) {
    stop("No DEG results found in the input data")
  }

  # Validate gene_sets input
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("gene_sets must be a named list where each element contains genes in a pathway")
  }

  # Initialize results storage
  pathway_results <- list()

  # Process each comparison
  for (comparison_name in names(deg_list)) {
    deg_data <- deg_list[[comparison_name]]

    # Filter DEGs by direction
    if (direction == "up") {
      genes <- deg_data$gene[deg_data$significant & deg_data$avg_log2FC > 0]
    } else if (direction == "down") {
      genes <- deg_data$gene[deg_data$significant & deg_data$avg_log2FC < 0]
    } else {
      genes <- deg_data$gene[deg_data$significant]
    }

    # Skip if no genes
    if (length(genes) == 0) {
      message(paste0("Skipping pathway analysis for ", comparison_name,
                     ": no ", direction, " regulated genes found"))
      next
    }

    # Run enrichment analysis
    enrichment_results <- simple_ora(
      gene_list = genes,
      gene_sets = gene_sets,
      universe = background,
      min_genes = min_genes,
      p_val_cutoff = p_adj_cutoff
    )

    # Store results
    if (nrow(enrichment_results) > 0) {
      pathway_results[[comparison_name]] <- list(
        enrichment = enrichment_results,
        genes = genes,
        direction = direction
      )
    }
  }

  # Return results
  return(list(
    pathway_results = pathway_results,
    metadata = list(
      min_genes = min_genes,
      p_adj_cutoff = p_adj_cutoff,
      direction = direction
    )
  ))
}
