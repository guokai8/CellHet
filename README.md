# CellHet: Advanced Cellular Heterogeneity Analysis for Single-Cell Data

[![GitHub issues](https://img.shields.io/github/issues/guokai8/CellHet)](https://github.com/guokai8/CellHet/issues)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

## Overview

CellHet (Cellular Heterogeneity) is a comprehensive R package designed to analyze differential gene expression and cellular heterogeneity in single-cell data. The package provides an integrated workflow for comparing cell types across conditions, identifying conserved and unique gene expression patterns, performing pathway enrichment analysis, integrating multi-modal data, and visualizing results through a rich set of plotting functions.

**Key Features:**
- Differential expression analysis across cell types and conditions
- Versatile visualizations with customizable heatmaps, barplots, dotplots, and UpSet plots
- Pathway enrichment analysis with visualization tools
- Multi-modal integration capabilities (RNA, ATAC, protein)
- Trajectory analysis and visualization
- Cell type response comparison and heterogeneity mapping
- Full compatibility with both Seurat and SingleCellExperiment objects

## Installation

```r
# Install from GitHub
devtools::install_github("guokai8/CellHet")

# Load the package
library(CellHet)
```

## Workflow Overview

![CellHet Workflow](workflow_diagram.png)

## Core Functionality

### Differential Expression Analysis

CellHet offers multiple approaches for DEG analysis with varying levels of complexity:

#### 1. Basic DEG Analysis

```r
# Compare DEGs across all cell types and conditions
deg_results <- compareDEGs(
  sce = seurat_object,              # Seurat or SingleCellExperiment object
  group_var = "condition",          # Group/condition column
  cell_type_var = "cell_type",      # Cell type column
  min_cells = 3,                    # Minimum cells per group
  test.use = "wilcox",              # Statistical test: wilcox, t, negbinom, etc.
  logfc.threshold = 0.25,           # Log fold-change threshold
  p_val_adj_threshold = 0.05,       # Adjusted p-value threshold
  custom_comparisons = NULL,        # Optional list of custom comparisons
  cores = 4                         # Number of cores for parallel processing
)
```

#### 2. Comprehensive Analysis

```r
# Complete analysis pipeline
full_results <- findDifferentialGenes(
  object = seurat_object,
  group_var = "condition",
  cell_type_var = "cell_type",
  min_cells = 3,
  logfc.threshold = 0.25,
  p_val_adj_threshold = 0.05,
  test.use = "wilcox",
  workers = 4,                       # Parallel processing
  gene_id_mapping = id_map_df,       # Optional gene ID mapping
  enrich = TRUE,                     # Run pathway analysis
  annot_data = annot_data,           # Annotation database
  use_gsea = FALSE,                  # ORA by default, set TRUE for GSEA
  return_summary = TRUE,             # Return summary visualization
  reference_group = "control"        # Use control as reference
)
```

#### 3. Targeted Comparison

```r
# Quick comparison of specific groups and cell types
t_cell_results <- quickCompare(
  object = seurat_object,
  cell_types = c("CD4+ T", "CD8+ T"),  # Specific cell types
  group1 = "treatment",                # First group
  group2 = "control",                  # Second group
  logfc.threshold = 0.25,
  p_val_adj_threshold = 0.05,
  return_plot = TRUE                   # Returns volcano plots
)
```

### Visualization

CellHet provides multiple visualization options to explore differential expression results:

#### Heatmaps

```r
# Create DEG heatmap
heatmap <- visualizeDEGHeatmap(
  deg_results = deg_results,
  comparison_col = "comparison",
  cell_type_col = "cell_type",
  up_color = "darkorange",
  down_color = "deepskyblue",
  show_counts = TRUE,                # Show count numbers on tiles
  group_regulation = TRUE,           # Group up/down regulation on x-axis
  show_reg_indicators = TRUE,        # Add color indicators for regulation
  cluster_rows = TRUE,               # Cluster cell types by similarity
  title = "DEGs Across Cell Types"
)
```

#### Barplots

```r
# Create barplots
barplot <- visualizeDEGBarplot(
  deg_results = deg_results,
  comparison_col = "comparison",
  cell_type_col = "cell_type",
  facet_by = "comparison",           # Create facets by comparison
  stacked = FALSE,                   # Side-by-side bars (TRUE for stacked)
  horizontal = FALSE,                # Vertical bars
  show_counts = TRUE                 # Show count numbers on bars
)
```

#### Dotplots

```r
# Create dotplot
dotplot <- visualizeDEGDotplot(
  deg_results = deg_results,
  comparison_col = "comparison",
  cell_type_col = "cell_type",
  count_display = "inside",          # Where to display counts: "none", "inside", "below"
  use_triangles = TRUE,              # Use triangles for up/down regulation
  offset = 0.2,                      # Horizontal offset between up/down points
  title = "DEG Dotplot"
)
```

#### UpSet Plots

```r
# Create UpSet plot for set intersection visualization
upset_plot <- visualizeDEGUpset(
  deg_results = deg_results,
  by_cell_type = TRUE,               # Show overlap across cell types
  direction = "up",                  # "up", "down", or "both"
  min_size = 5,                      # Minimum intersection size
  highlight_exclusive = TRUE,        # Highlight cell type-specific DEGs
  highlight_shared = TRUE,           # Highlight DEGs shared by all
  title = "Shared and Unique Up-regulated Genes"
)
```

#### Comprehensive Visualization

```r
# Create multi-panel visualization
pattern_viz <- visualizeDEGPatterns(
  deg_results = deg_results,
  plot_types = c("heatmap", "barplot", "dotplot"),
  ncol = 1,                          # Number of columns for layout
  title = "DEG Patterns Across Cell Types",
  cluster_rows = TRUE,
  hide_empty_cell_types = TRUE       # Hide cell types with no DEGs
)
```

### Pathway Analysis

Perform pathway enrichment analysis on DEG results:

```r
# Run pathway analysis
pathway_results <- runPathwayAnalysis(
  deg_results = deg_results,
  annot_data = go_annot,             # Annotation data from richR package
  annot_type = "GO",                 # "GO", "KEGG", "KEGGM", or "MSIGDB"
  ontology = "BP",                   # For GO: "BP", "MF", or "CC"
  direction = "both",                # Consider up/down/both DEGs
  min_genes = 5,                     # Min genes per pathway
  use_gsea = FALSE                   # Use ORA (FALSE) or GSEA (TRUE)
)

# Visualize pathway results
pathway_plot <- visualizePathwaysDotplot(
  pathway_results = pathway_results,
  top_n = 10,                        # Top pathways to show
  color_by = "RichFactor",           # Color by: "RichFactor", "Padj", etc.
  size_by = "Pvalue",                # Size by: "Pvalue", "Padj", etc.
  flip = TRUE,                       # Flip coordinates
  title = "Pathway Enrichment"
)
```

### Multi-Modal Integration

Integrate and analyze data across multiple modalities:

```r
# Integrate data across modalities
integration_results <- integrateMultiModalData(
  data_list = list(rna_data, atac_data, protein_data),
  modality_names = c("RNA", "ATAC", "Protein"),
  cell_type_var = "cell_type",
  integration_method = "harmony",    # "harmony", "canonical_correlation", "mnn", "seurat"
  dim_reduction = "umap",            # Dimension reduction method
  n_dims = 30,                       # Dimensions for integration
  n_features = 2000                  # Features per modality
)

# Visualize integration results
integration_viz <- visualizeIntegration(
  integration_results = integration_results,
  color_by = c("modality", "cell_type"),
  features = top_variable_features,  # Features to visualize
  plot_types = c("dimplot", "feature", "cluster", "metrics", "density"),
  split_by = "condition",
  interactive = FALSE
)

# Optimize integration parameters
optimized_integration <- optimizeIntegrationParameters(
  data_list = list(rna_data, atac_data),
  modality_names = c("RNA", "ATAC"),
  cell_type_var = "cell_type",
  optimization_metric = "mixing",    # "mixing", "silhouette", "kbet", "lisi"
  methods = c("harmony", "mnn"),     # Methods to test
  n_dims_range = c(10, 20, 30),      # Dimensions to test
  n_features_range = c(1000, 2000, 3000)
)
```

### Advanced Analyses

#### Trajectory Analysis

```r
# Analyze gene expression along a trajectory
trajectory_results <- analyzeDEGTrajectory(
  sce = seurat_object,
  trajectory_var = "pseudotime",     # Column with trajectory values
  cell_type_var = "cell_type",
  group_var = "condition",           # Optional grouping
  n_bins = 10,                       # Number of trajectory bins
  smooth_method = "loess",           # "loess", "gam", or "none"
  n_top_genes = 50,                  # Top dynamic genes
  cluster_genes = TRUE               # Cluster genes by pattern
)

# Plot trajectories for specific genes
gene_trajectory <- plotGeneTrajectory(
  trajectory_results = trajectory_results,
  genes = c("CD3D", "MS4A1", "NKG7"),
  plot_type = "line",                # "line", "box", or "violin"
  add_loess_fit = TRUE,
  facet_by_cluster = TRUE
)

# Compare trajectories between conditions
diff_trajectory <- diffTrajectoryAnalysis(
  sce = seurat_object,
  trajectory_var = "pseudotime",
  group_var = "condition",
  cell_type_var = "cell_type",
  n_bins = 10,
  test_method = "anova",             # Statistical method
  smooth_method = "loess"
)
```

#### Cell Type Response Comparison

```r
# Compare cell type responses to a reference
response_comparison <- compareResponses(
  deg_results = deg_results,
  reference_cell_type = "CD4+ T",    # Reference cell type
  plot_type = "heatmap",             # "heatmap", "barplot", "radar"
  similarity_metric = "jaccard"      # Similarity metric
)

# Visualize responses across conditions
response_viz <- visualizeReferenceComparisons(
  deg_results = deg_results,
  reference_group = "control",       # Reference condition
  plot_type = "all",                 # "heatmap", "similarity", "trajectory", "all"
  interactive = FALSE
)
```

#### Heterogeneity Mapping

```r
# Create heterogeneity map
het_map <- visualizeHeterogeneity(
  deg_results = deg_results,
  pathway_results = pathway_results, # Optional
  layout_type = "network",           # "heatmap", "network", "umap", "circular"
  level = "both",                    # "gene", "pathway", "both"
  similarity_metric = "jaccard"
)
```

### Data Export

Export results for sharing and further analysis:

```r
# Export to Excel
exportDEGResults(
  deg_results = deg_results,
  file_path = "deg_results.xlsx",
  split_by = "cell_type",           # How to organize sheets
  include_pathways = TRUE           # Include pathway results
)
```

## Example Workflow with Simulated Data

CellHet includes functions to generate simulated data, allowing users to test the functionality without their own datasets. The following example demonstrates a complete workflow using simulated data:

```r
# Load required libraries
library(CellHet)
library(Seurat)
library(dplyr)

# Generate simulated single-cell data
set.seed(42) # For reproducibility
sim_data <- simulateSingleCellData(
  n_cells = 1000,              # 1000 cells
  n_genes = 2000,              # 2000 genes
  n_cell_types = 5,            # 5 cell types
  n_conditions = 2,            # 2 conditions (e.g., treatment vs control)
  pct_de_genes = 0.1,          # 10% of genes are differentially expressed
  effect_size = 1.5,           # Effect size for DE genes
  batch_effect = TRUE,         # Include batch effects
  add_trajectory = TRUE,       # Add pseudotime trajectory
  return_type = "seurat"       # Return a Seurat object
)

# Explore the simulated data
print(table(sim_data$cell_type, sim_data$condition))
print(head(sim_data@meta.data))

# Run differential expression analysis
deg_results <- findDifferentialGenes(
  object = sim_data,
  group_var = "condition",
  cell_type_var = "cell_type",
  reference_group = "Condition1",
  workers = 4,
  return_summary = TRUE
)

# Visualize results
heatmap <- visualizeDEGHeatmap(deg_results)
barplot <- visualizeDEGBarplot(deg_results)
upset <- visualizeDEGUpset(deg_results, by_cell_type = TRUE)

# Compare responses across cell types
response_viz <- visualizeReferenceComparisons(
  deg_results = deg_results,
  reference_group = "Condition1",
  plot_type = "all"
)

# Generate simulated multi-modal data
multi_modal_data <- simulateMultiModalData(
  n_cells = 800,
  n_modalities = 2,                          # RNA and ATAC
  n_features = c(2000, 1000),                # 2000 genes, 1000 peaks
  n_cell_types = 4,                          # 4 cell types
  n_conditions = 2,                          # 2 conditions
  shared_cell_types = TRUE,                  # Same cell types across modalities
  cell_overlap = 0.7,                        # 70% of cells shared between modalities
  feature_overlap = 0.3,                     # 30% of features have correlated expression
  correlation_strength = 0.7,                # Correlation strength between modalities
  return_type = "seurat"                     # Return Seurat objects
)

# Access individual modalities
rna_data <- multi_modal_data$Modality1
atac_data <- multi_modal_data$Modality2

# Integrate the modalities
integration_results <- integrateMultiModalData(
  data_list = multi_modal_data,
  modality_names = c("RNA", "ATAC"),
  cell_type_var = "cell_type",
  integration_method = "harmony",
  dim_reduction = "umap",
  n_dims = 30
)

# Visualize integration
integration_viz <- visualizeIntegration(
  integration_results = integration_results,
  color_by = c("modality", "cell_type", "condition"),
  plot_types = c("dimplot", "density")
)

# Analyze gene expression trajectories using pseudotime
trajectory_results <- analyzeDEGTrajectory(
  sce = sim_data,
  trajectory_var = "pseudotime",
  cell_type_var = "cell_type",
  n_bins = 10,
  smooth_method = "loess"
)

# Export results
temp_file <- tempfile(fileext = ".xlsx")
exportDEGResults(
  deg_results = deg_results,
  file_path = temp_file,
  split_by = "cell_type"
)
```

## Working with Simulated Data

CellHet provides functions to generate simulated data for testing and method development:

```r
# Simulate single-cell data
sim_data <- simulateSingleCellData(
  n_cells = 1000,
  n_genes = 2000,
  n_cell_types = 5,
  n_conditions = 2,
  return_type = "seurat",
  batch_effect = TRUE,
  add_trajectory = TRUE
)

# Simulate multi-modal data
multi_modal_data <- simulateMultiModalData(
  n_cells = 1000,
  n_modalities = 3,
  n_features = c(2000, 1500, 500),
  n_cell_types = 5,
  n_conditions = 2,
  return_type = "seurat",
  shared_cell_types = TRUE,
  cell_overlap = 0.7,
  feature_overlap = 0.3
)
```

## Citation

If you use CellHet in your research, please cite:

[(https://github.com/guokai8/CellHet)]

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

For questions, issues, or feature requests, please open an issue on the GitHub repository.
