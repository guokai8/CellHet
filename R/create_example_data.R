#' Create Example Single-Cell Dataset
#'
#' Creates a small example Seurat object for testing and demonstration purposes.
#' The dataset contains simulated PBMC expression data with realistic gene names,
#' cell type annotations, conditions, and pseudotime values.
#'
#' @param n_cells Number of cells to simulate (default: 500)
#' @param n_genes Number of genes to simulate (default: 200)
#' @param n_cell_types Number of cell types (default: 3, max: 6)
#' @param n_conditions Number of conditions (default: 3)
#' @param seed Random seed for reproducibility (default: 123)
#'
#' @return A Seurat object with:
#'   \itemize{
#'     \item Normalized expression data with realistic PBMC gene names
#'     \item Cell type annotations (cell_type): CD4 T cells, CD8 T cells, B cells, NK cells, Monocytes, etc.
#'     \item Condition labels (condition)
#'     \item Pseudotime values (pseudotime)
#'     \item Cluster assignments (seurat_clusters)
#'   }
#'
#' @examples
#' # Create example dataset
#' pbmc_example <- createExampleData()
#'
#' # Check structure
#' pbmc_example
#'
#' # View metadata
#' head(pbmc_example@meta.data)
#'
#' @export
createExampleData <- function(n_cells = 500,
                              n_genes = 200,
                              n_cell_types = 3,
                              n_conditions = 3,
                              seed = 123) {

  # Check if Seurat is available
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required to create example data. Please install it with: install.packages('Seurat')")
  }

  set.seed(seed)

  # Realistic PBMC gene names (common markers and housekeeping genes)
  pbmc_genes <- c(
    # T cell markers
    "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "IL7R", "CCR7", "SELL", "LEF1",
    "CD28", "LCK", "ZAP70", "LAT", "GZMK", "GZMB", "PRF1", "GNLY", "NKG7", "CCL5",

    # B cell markers
    "CD19", "MS4A1", "CD79A", "CD79B", "IGHM", "IGHD", "IGHA1", "IGHA2", "IGHG1", "IGHG2",
    "CD27", "TCL1A", "FCER2", "PAX5", "CD38", "SDC1", "XBP1", "PRDM1", "IRF4", "CD22",

    # NK cell markers
    "NCAM1", "FCGR3A", "KLRB1", "KLRD1", "KLRF1", "KLRC1", "KLRC2", "NCR1", "NCR3", "TYROBP",

    # Monocyte/Macrophage markers
    "CD14", "CD16", "FCGR1A", "CSF1R", "CD68", "LYZ", "S100A8", "S100A9", "S100A12", "VCAN",
    "CD163", "MSR1", "MRC1", "CD86", "IL1B", "TNF", "CCL2", "CCL3", "CCL4", "CXCL8",

    # Dendritic cell markers
    "CD1C", "FCER1A", "CLEC10A", "CLEC9A", "XCR1", "BATF3", "IRF8", "FLT3", "CD209", "CLEC4C",

    # Housekeeping and common genes
    "ACTB", "GAPDH", "B2M", "RPL13", "RPS18", "UBC", "PPIA", "HPRT1", "TBP", "YWHAZ",
    "RPLP0", "RPS27A", "RPL32", "RPS29", "TUBB", "EEF1A1", "HSP90AB1", "HSPA8", "VIM", "TPT1",

    # Cytokines and signaling
    "IFNG", "IL2", "IL4", "IL6", "IL10", "IL12A", "IL12B", "IL17A", "IL21", "IL23A",
    "TGFB1", "CSF2", "CSF3", "FAS", "FASLG", "TNFSF10", "CD40LG", "CD40", "CTLA4", "PDCD1",

    # Transcription factors
    "FOXP3", "TBX21", "GATA3", "RORC", "BCL6", "STAT1", "STAT3", "STAT4", "STAT5A", "STAT6",
    "NFKB1", "NFKB2", "REL", "RELA", "RELB", "JUN", "FOS", "MYC", "TP53", "RUNX1",

    # Proliferation markers
    "MKI67", "PCNA", "TOP2A", "CDK1", "CCNA2", "CCNB1", "CCND1", "CCNE1", "CDK2", "CDK4",

    # Adhesion and migration
    "ITGAL", "ITGAM", "ITGAX", "ITGB2", "SELP", "SELE", "PECAM1", "VCAM1", "ICAM1", "CXCR4",
    "CXCR3", "CCR2", "CCR4", "CCR5", "CCR6", "CCR9", "CXCL9", "CXCL10", "CXCL11", "CXCL12"
  )

  # Extend gene list if needed
  if (n_genes > length(pbmc_genes)) {
    additional_genes <- paste0("GENE", 1:(n_genes - length(pbmc_genes)))
    pbmc_genes <- c(pbmc_genes, additional_genes)
  }

  # Select genes to use
  selected_genes <- pbmc_genes[1:n_genes]

  # Create simulated expression matrix
  counts <- matrix(
    rpois(n_cells * n_genes, lambda = 5),
    nrow = n_genes,
    ncol = n_cells
  )

  # Add some variation
  counts <- counts + matrix(
    rnorm(n_cells * n_genes, mean = 0, sd = 2),
    nrow = n_genes,
    ncol = n_cells
  )
  counts[counts < 0] <- 0

  # Use realistic gene names
  rownames(counts) <- selected_genes
  colnames(counts) <- paste0("Cell-", 1:n_cells)

  # Realistic PBMC cell types
  pbmc_cell_types <- c("CD4_T", "CD8_T", "B_cells", "NK_cells", "Monocytes", "DCs")
  selected_cell_types <- pbmc_cell_types[1:min(n_cell_types, 6)]

  # Create Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts,
    project = "CellHet_PBMC_Example"
  )

  # Add metadata with realistic cell types
  seurat_obj$cell_type <- sample(
    selected_cell_types,
    n_cells,
    replace = TRUE
  )

  seurat_obj$condition <- sample(
    c("Control", "Stimulated", "Treatment")[1:n_conditions],
    n_cells,
    replace = TRUE
  )

  seurat_obj$pseudotime <- runif(n_cells, min = 0, max = 1)

  seurat_obj$seurat_clusters <- as.factor(
    sample(0:(n_cell_types - 1), n_cells, replace = TRUE)
  )

  # Normalize data
  seurat_obj <- Seurat::NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    verbose = FALSE
  )

  # Add cell-type specific marker expression
  # Define marker genes for each cell type
  marker_genes <- list(
    CD4_T = c("CD3D", "CD3E", "CD4", "IL7R", "CCR7", "LEF1", "SELL", "CD28", "LCK", "TCF7"),
    CD8_T = c("CD3D", "CD3E", "CD8A", "CD8B", "GZMK", "GZMB", "PRF1", "CCL5", "NKG7", "GNLY"),
    B_cells = c("CD19", "MS4A1", "CD79A", "CD79B", "IGHM", "IGHD", "CD27", "TCL1A", "PAX5", "CD22"),
    NK_cells = c("NCAM1", "FCGR3A", "KLRB1", "KLRD1", "KLRF1", "NKG7", "GNLY", "PRF1", "GZMB", "TYROBP"),
    Monocytes = c("CD14", "CD16", "FCGR1A", "LYZ", "S100A8", "S100A9", "CD68", "CSF1R", "VCAN", "CCL2"),
    DCs = c("CD1C", "FCER1A", "CLEC10A", "CLEC9A", "FLT3", "IRF8", "BATF3", "CD209", "XCR1", "CLEC4C")
  )

  # Boost expression of marker genes for each cell type
  for (ct_idx in 1:length(selected_cell_types)) {
    ct_name <- selected_cell_types[ct_idx]
    cell_idx <- which(seurat_obj$cell_type == ct_name)

    if (length(cell_idx) > 0 && ct_name %in% names(marker_genes)) {
      # Get marker genes for this cell type that are in our gene set
      ct_markers <- marker_genes[[ct_name]]
      ct_markers <- ct_markers[ct_markers %in% selected_genes]

      if (length(ct_markers) > 0) {
        marker_idx <- which(rownames(counts) %in% ct_markers)
        if (length(marker_idx) > 0) {
          # Boost marker gene expression 3-5x
          counts[marker_idx, cell_idx] <- counts[marker_idx, cell_idx] * runif(length(marker_idx), 3, 5)
        }
      }
    }
  }

  # Add condition-specific differential expression
  # Upregulate some genes in "Stimulated" or "Treatment" conditions
  stim_genes <- c("IFNG", "TNF", "IL2", "IL6", "FOS", "JUN", "NFKB1", "STAT1", "STAT3", "MYC")
  stim_genes <- stim_genes[stim_genes %in% selected_genes]

  if (length(stim_genes) > 0 && "Stimulated" %in% unique(seurat_obj$condition)) {
    stim_cells <- which(seurat_obj$condition == "Stimulated")
    if (length(stim_cells) > 0) {
      stim_gene_idx <- which(rownames(counts) %in% stim_genes)
      if (length(stim_gene_idx) > 0) {
        counts[stim_gene_idx, stim_cells] <- counts[stim_gene_idx, stim_cells] * runif(length(stim_gene_idx), 2, 4)
      }
    }
  }

  # Update with modified counts
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts,
    project = "CellHet_PBMC_Example",
    meta.data = seurat_obj@meta.data
  )

  seurat_obj <- Seurat::NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    verbose = FALSE
  )

  message("Created example PBMC Seurat object with:")
  message("  - ", n_cells, " cells")
  message("  - ", n_genes, " genes (realistic PBMC gene names)")
  message("  - ", n_cell_types, " cell types: ", paste(selected_cell_types, collapse = ", "))
  message("  - ", n_conditions, " conditions")

  return(seurat_obj)
}
