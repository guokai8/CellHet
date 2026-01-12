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

## Quick Start

Here's a complete, runnable example using realistic PBMC data with human gene symbols:

```r
library(CellHet)
library(Seurat)

# Step 1: Create example PBMC data with real human gene symbols
# This creates a Seurat object with realistic markers like CD3D, CD4, CD8A, etc.
pbmc <- createExampleData(
  n_cells = 500,
  n_genes = 200,
  n_cell_types = 4,
  n_conditions = 2,
  seed = 123
)

# View the data structure
print(pbmc)
head(pbmc@meta.data)
table(pbmc$cell_type, pbmc$condition)

# Step 2: Run differential expression analysis
deg_results <- findDifferentialGenes(
  object = pbmc,
  group_var = "condition",
  cell_type_var = "cell_type",
  reference_group = "Control",
  logfc_threshold = 0.25,
  p_val_adj_threshold = 0.05,
  workers = 2
)

# Step 3: View results
# Check summary
print(deg_results$summary)

# View top DEGs for a specific cell type and comparison
head(deg_results$all_degs)

# Step 4: Visualize results
# Heatmap showing DEG counts
heatmap_plot <- visualizeDEGHeatmap(
  deg_results = deg_results,
  show_counts = TRUE,
  cluster_rows = TRUE
)
print(heatmap_plot)

# Barplot
barplot <- visualizeDEGBarplot(
  deg_results = deg_results,
  facet_by = "comparison"
)
print(barplot)

# UpSet plot showing shared genes
upset_plot <- visualizeDEGUpset(
  deg_results = deg_results,
  by_cell_type = TRUE,
  direction = "both"
)
print(upset_plot)

# Step 5: Export results
exportDEGResults(
  deg_results = deg_results,
  file_path = "pbmc_deg_results.xlsx",
  split_by = "cell_type"
)
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
  object = pbmc,                    # Seurat or SingleCellExperiment object
  group_var = "condition",          # Group/condition column
  cell_type_var = "cell_type",      # Cell type column
  min_cells_per_group = 3,          # Minimum cells per group
  test_method = "wilcox",           # Statistical test: wilcox, t, negbinom, etc.
  logfc_threshold = 0.25,           # Log fold-change threshold
  p_val_adj_threshold = 0.05,       # Adjusted p-value threshold
  custom_comparisons = NULL,        # Optional list of custom comparisons
  cores = 4                         # Number of cores for parallel processing
)
```

#### 2. Comprehensive Analysis

```r
# Complete analysis pipeline with pathway enrichment
full_results <- findDifferentialGenes(
  object = pbmc,
  group_var = "condition",
  cell_type_var = "cell_type",
  min_cells_per_group = 3,
  logfc_threshold = 0.25,
  p_val_adj_threshold = 0.05,
  test_method = "wilcox",
  workers = 2,                       # Parallel processing
  gene_id_mapping = NULL,            # Optional gene ID mapping
  enrich = FALSE,                    # Set TRUE to run pathway analysis
  annot_data = NULL,                 # Annotation database (e.g., from richR)
  use_gsea = FALSE,                  # ORA by default, set TRUE for GSEA
  return_summary = TRUE,             # Return summary visualization
  reference_group = "Control"        # Use Control as reference
)
```

#### 3. Targeted Comparison

```r
# Quick comparison of specific groups and cell types
t_cell_results <- quickCompare(
  object = pbmc,
  cell_types = c("CD4_T", "CD8_T"),    # Specific cell types
  group1 = "Stimulated",               # First group
  group2 = "Control",                  # Second group
  logfc_threshold = 0.25,
  p_val_adj_threshold = 0.05,
  test_method = "wilcox",
  return_plot = TRUE                   # Returns volcano plots
)

# Access results
print(t_cell_results$results)
print(t_cell_results$plots$CD4_T)    # View CD4 T cell volcano plot
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
# Note: Requires richR package and annotation database
# Example with GO enrichment
\dontrun{
  # First get annotation data
  library(richR)
  # For human genes
  go_annot <- getGO(organism = "human")

  # Run pathway analysis on DEG results
  pathway_results <- runPathwayAnalysis(
    deg_results = deg_results,
    annot_data = go_annot,           # Annotation data from richR package
    annot_type = "GO",               # "GO", "KEGG", "KEGGM", or "MSIGDB"
    ontology = "BP",                 # For GO: "BP", "MF", or "CC"
    direction = "both",              # Consider up/down/both DEGs
    min_genes = 5,                   # Min genes per pathway
    use_gsea = FALSE                 # Use ORA (FALSE) or GSEA (TRUE)
  )

  # Visualize pathway results
  pathway_plot <- visualizePathwaysDotplot(
    pathway_results = pathway_results,
    top_n = 10,                      # Top pathways to show
    color_by = "RichFactor",         # Color by: "RichFactor", "Padj", etc.
    size_by = "Pvalue",              # Size by: "Pvalue", "Padj", etc.
    flip = TRUE,                     # Flip coordinates
    title = "Pathway Enrichment"
  )
}
```

### Multi-Modal Integration

Integrate and analyze data across multiple modalities:

```r
# First, create multi-modal data
multi_modal_data <- simulateMultiModalData(
  n_cells = 800,
  n_modalities = 2,
  n_features = c(2000, 1000),
  n_cell_types = 4,
  n_conditions = 2,
  return_type = "seurat"
)

# Integrate data across modalities
integration_results <- integrateMultiModalData(
  data_list = multi_modal_data,
  modality_names = c("RNA", "ATAC"),
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
  plot_types = c("dimplot", "density"),
  interactive = FALSE
)

# Optimize integration parameters (test different settings)
\dontrun{
  optimized_integration <- optimizeIntegrationParameters(
    data_list = multi_modal_data,
    modality_names = c("RNA", "ATAC"),
    cell_type_var = "cell_type",
    optimization_metric = "mixing",    # "mixing", "silhouette", "kbet", "lisi"
    methods = c("harmony", "mnn"),     # Methods to test
    n_dims_range = c(10, 20, 30),      # Dimensions to test
    n_features_range = c(1000, 2000, 3000)
  )
}
```

### Advanced Analyses

#### Trajectory Analysis

```r
# Note: Pseudotime must be calculated first (requires monocle3)
# See the Pseudotime Analysis section for calculatePseudotime()

# Analyze gene expression along a trajectory
\dontrun{
  trajectory_results <- analyzeDEGTrajectory(
    sce = pbmc,
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
    sce = pbmc,
    trajectory_var = "pseudotime",
    group_var = "condition",
    cell_type_var = "cell_type",
    n_bins = 10,
    test_method = "anova",             # Statistical method
    smooth_method = "loess"
  )
}
```

#### Cell Type Response Comparison

```r
# Compare cell type responses to a reference
response_comparison <- compareResponses(
  deg_results = deg_results,
  reference_cell_type = "CD4_T",     # Reference cell type
  plot_type = "heatmap",             # "heatmap", "barplot", "radar"
  similarity_metric = "jaccard"      # Similarity metric
)

# Visualize responses across conditions
response_viz <- visualizeReferenceComparisons(
  deg_results = deg_results,
  reference_group = "Control",       # Reference condition
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

## Complete Workflow Examples

### Example 1: Comprehensive DEG Analysis with Real Gene Symbols

Building on the Quick Start example, here's a more comprehensive analysis:

```r
library(CellHet)
library(Seurat)
library(dplyr)

# Create example PBMC data with realistic human gene symbols
pbmc <- createExampleData(
  n_cells = 500,
  n_genes = 200,
  n_cell_types = 4,
  n_conditions = 2,
  seed = 123
)

# View gene names (real human gene symbols like CD3D, CD4, CD8A, etc.)
head(rownames(pbmc))

# Check the data
print(table(pbmc$cell_type, pbmc$condition))

# Run comprehensive differential expression analysis
deg_results <- findDifferentialGenes(
  object = pbmc,
  group_var = "condition",
  cell_type_var = "cell_type",
  reference_group = "Control",
  logfc_threshold = 0.25,
  p_val_adj_threshold = 0.05,
  workers = 2,
  return_summary = TRUE
)

# View summary
print(deg_results$summary)

# View specific DEGs
cd4_degs <- deg_results$all_degs %>%
  filter(cell_type == "CD4_T", significant == TRUE) %>%
  arrange(desc(abs(avg_log2FC)))
print(head(cd4_degs))

# Create comprehensive visualizations
heatmap <- visualizeDEGHeatmap(
  deg_results = deg_results,
  show_counts = TRUE,
  cluster_rows = TRUE,
  title = "DEGs Across PBMC Cell Types"
)

barplot <- visualizeDEGBarplot(
  deg_results = deg_results,
  facet_by = "comparison",
  show_counts = TRUE
)

dotplot <- visualizeDEGDotplot(
  deg_results = deg_results,
  use_triangles = TRUE,
  count_display = "inside"
)

upset_plot <- visualizeDEGUpset(
  deg_results = deg_results,
  by_cell_type = TRUE,
  direction = "both",
  highlight_exclusive = TRUE
)

# Analyze shared DEGs
shared_analysis <- findSharedDEGs(
  deg_results = deg_results,
  min_deg_count = 5,
  plot_type = "upset",
  direction = "both"
)

# View shared genes
print(shared_analysis$gene_sets)

# Create multi-panel visualization
pattern_viz <- visualizeDEGPatterns(
  deg_results = deg_results,
  plot_types = c("heatmap", "barplot", "dotplot"),
  ncol = 1
)

# Export to Excel
exportDEGResults(
  deg_results = deg_results,
  file_path = "pbmc_deg_results.xlsx",
  split_by = "cell_type"
)
```

### Example 2: Large-Scale Simulated Study

For testing with larger datasets or custom parameters:

```r
# Generate larger simulated dataset
set.seed(42)
sim_data <- simulateSingleCellData(
  n_cells = 1000,              # 1000 cells
  n_genes = 2000,              # 2000 genes
  n_cell_types = 5,            # 5 cell types
  n_conditions = 2,            # 2 conditions
  pct_de_genes = 0.1,          # 10% DE genes
  effect_size = 1.5,           # Effect size
  batch_effect = TRUE,         # Include batch effects
  add_trajectory = TRUE,       # Add pseudotime
  return_type = "seurat"
)

# Run analysis
deg_results <- findDifferentialGenes(
  object = sim_data,
  group_var = "condition",
  cell_type_var = "cell_type",
  reference_group = "Condition1",
  workers = 4,
  return_summary = TRUE
)

# Visualize results
visualizeDEGPatterns(deg_results, plot_types = c("heatmap", "barplot"))
```

### Example 3: Multi-Modal Integration

```r
# Generate multi-modal data (RNA + ATAC)
multi_modal_data <- simulateMultiModalData(
  n_cells = 800,
  n_modalities = 2,
  n_features = c(2000, 1000),
  n_cell_types = 4,
  n_conditions = 2,
  shared_cell_types = TRUE,
  cell_overlap = 0.7,
  feature_overlap = 0.3,
  correlation_strength = 0.7,
  return_type = "seurat"
)

# Integrate modalities
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
  color_by = c("modality", "cell_type"),
  plot_types = c("dimplot", "density")
)

# Evaluate integration quality
eval_results <- evaluateIntegration(
  integration_results = integration_results,
  cell_type_var = "cell_type",
  metrics = c("silhouette", "mixing")
)
print(eval_results$metrics)
```

### Example 4: Pseudotime Analysis

```r
# Using the pbmc data from Quick Start
pbmc <- createExampleData(n_cells = 500, n_genes = 200, seed = 123)

# Note: For real pseudotime analysis, you would need monocle3 installed
# This example shows the function usage
\dontrun{
  # Calculate pseudotime
  pbmc <- calculatePseudotime(
    pbmc,
    root_type = "CD4_T",
    cell_type_var = "cell_type",
    reduction = "umap"
  )

  # Plot pseudotime
  plotPseudotime(pbmc, color_by = "pseudotime")
  plotPseudotime(pbmc, color_by = "cell_type")

  # Get statistics
  stats <- getPseudotimeStats(pbmc, group_by = "cell_type")
  print(stats)

  # Analyze trajectory
  trajectory_results <- analyzeDEGTrajectory(
    sce = pbmc,
    trajectory_var = "pseudotime",
    cell_type_var = "cell_type",
    n_bins = 10,
    smooth_method = "loess"
  )
}
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

# Create a quick example dataset with realistic PBMC markers
pbmc_example <- createExampleData(
  n_cells = 500,
  n_genes = 200,
  n_cell_types = 3,
  n_conditions = 3,
  seed = 123
)

# View metadata
head(pbmc_example@meta.data)
```

## Additional Analysis Functions

### Pseudotime Analysis

Calculate and visualize pseudotime trajectories using Monocle3:

```r
# Calculate pseudotime
seurat_obj <- calculatePseudotime(
  seurat_obj,
  root_type = "Stem_Cells",
  cell_type_var = "cell_type",
  reduction = "umap"
)

# Plot pseudotime trajectory
plotPseudotime(seurat_obj, color_by = "pseudotime")
plotPseudotime(seurat_obj, color_by = "cell_type")

# Get pseudotime statistics
getPseudotimeStats(seurat_obj)
getPseudotimeStats(seurat_obj, group_by = "cell_type")
```

### DEG Analysis Utilities

#### Summarize DEG Results

```r
# Create summary visualizations
summary_barplot <- summarizeDEGs(
  deg_results = deg_results,
  plot_type = "barplot",
  direction = "both",
  interactive = FALSE
)

summary_heatmap <- summarizeDEGs(
  deg_results = deg_results,
  plot_type = "heatmap",
  direction = "up"
)

summary_dotplot <- summarizeDEGs(
  deg_results = deg_results,
  plot_type = "dotplot",
  direction = "down"
)
```

#### Find Shared DEGs

```r
# Analyze shared and unique DEGs across comparisons
shared_degs <- findSharedDEGs(
  deg_results = deg_results,
  cell_types = c("CD4+ T", "CD8+ T"),
  min_deg_count = 10,
  plot_type = "upset",
  direction = "both"
)

# Access shared gene sets
shared_degs$gene_sets
shared_degs$overlaps
shared_degs$plot

# Create different visualizations
venn_results <- findSharedDEGs(deg_results, plot_type = "venn")
heatmap_results <- findSharedDEGs(deg_results, plot_type = "heatmap")
network_results <- findSharedDEGs(deg_results, plot_type = "network")
```

#### Convert External DEG Results

```r
# Convert Seurat FindMarkers results to CellHet format
# Assuming you have DEG results from Seurat
deg_list <- list(
  "CD4_T_KO_vs_WT" = seurat_deg_results_cd4,
  "CD8_T_KO_vs_WT" = seurat_deg_results_cd8
)

cellhet_deg <- convertToCellHetDEG(
  deg_list = deg_list,
  p_val_adj_threshold = 0.05,
  logfc_threshold = 0.25
)

# Now use with CellHet visualization functions
visualizeDEGHeatmap(cellhet_deg)
```

### Multi-Modal Integration Functions

#### Find Conserved Features

```r
# Find features conserved across modalities
conserved_features <- findConservedFeatures(
  integration_results = integration_results,
  modality_names = c("RNA", "ATAC"),
  n_features = 100,
  correlation_threshold = 0.5
)
```

#### Find Multi-Modal Signatures

```r
# Identify cell type signatures across modalities
signatures <- findMultiModalSignatures(
  integration_results = integration_results,
  cell_type_var = "cell_type",
  modalities = c("RNA", "Protein"),
  n_features_per_modality = 50
)
```

#### Calculate Modality Agreement

```r
# Assess agreement between modalities
agreement <- calculateModalityAgreement(
  integration_results = integration_results,
  cell_type_var = "cell_type",
  metric = "correlation"
)
```

#### Evaluate Integration Quality

```r
# Evaluate integration performance
eval_results <- evaluateIntegration(
  integration_results = integration_results,
  cell_type_var = "cell_type",
  metrics = c("silhouette", "mixing", "kbet", "lisi")
)

# View evaluation metrics
print(eval_results$metrics)
plot(eval_results$plots$silhouette)
```

#### Feature Importance Analysis

```r
# Analyze feature importance in integration
feature_importance <- analyzeFeatureImportance(
  integration_results = integration_results,
  method = "random_forest",
  n_top_features = 50
)

# Plot top features
plot(feature_importance$plot)
```

#### Save and Load Integration Models

```r
# Export integration model for reuse
exportImportIntegrationModel(
  integration_results = integration_results,
  file_path = "integration_model.rds",
  mode = "export"
)

# Import saved model
loaded_model <- exportImportIntegrationModel(
  file_path = "integration_model.rds",
  mode = "import"
)

# Apply model to new data
projected_data <- applyIntegrationModel(
  new_data = new_seurat_obj,
  integration_model = loaded_model
)

# Project new data using existing integration
projected_results <- projectNewData(
  new_data = new_cells,
  integration_results = integration_results,
  reference_data = reference_data
)
```

### Enrichment Analysis Utilities

#### richR Integration

```r
# Run enrichment using richR package
richr_results <- runRichREnrichment(
  gene_list = deg_genes,
  organism = "human",
  database = "GO",
  ontology = "BP",
  pvalue_cutoff = 0.05
)
```

#### Export Enrichment Results

```r
# Export pathway results to Excel
exportEnrichmentResults(
  pathway_results = pathway_results,
  file_path = "pathway_results.xlsx",
  split_by = "cell_type"
)
```

#### Save Pathway Dotplots

```r
# Save pathway dotplot as high-resolution image
savePathwayDotplot(
  pathway_plot = pathway_plot,
  file_path = "pathway_dotplot.pdf",
  width = 10,
  height = 8,
  dpi = 300
)
```

### Color Palettes

CellHet provides a consistent color palette for all visualizations:

```r
# Get CellHet colors
colors <- cellhet_colors()         # All colors
colors_10 <- cellhet_colors(10)    # First 10 colors
colors_alpha <- cellhet_colors(5, alpha = 0.7)  # With transparency

# Use in ggplot2
library(ggplot2)
ggplot(data, aes(x = x, y = y, color = group)) +
  geom_point() +
  scale_color_cellhet()

ggplot(data, aes(x = x, y = y, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_cellhet()
```

### Interactive Shiny Application

Launch the interactive CellHet Shiny app for exploratory analysis:

```r
# Launch the app
runCellHetApp()

# Launch on specific port
runCellHetApp(port = 3838)

# Launch without browser auto-open
runCellHetApp(launch.browser = FALSE)
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
