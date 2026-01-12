# CellHet: Cellular Heterogeneity Analysis for Single-Cell Data

[![GitHub issues](https://img.shields.io/github/issues/guokai8/CellHet)](https://github.com/guokai8/CellHet/issues)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

## Overview

CellHet is an R package for analyzing cellular heterogeneity in single-cell data. It provides tools for differential expression analysis, pathway enrichment, multi-modal integration, and visualization.

**Key Features:**
- Differential expression analysis across cell types and conditions
- Rich visualization suite (heatmaps, barplots, dotplots, UpSet plots)
- Pathway enrichment analysis
- Multi-modal data integration (RNA, ATAC, protein)
- Pseudotime trajectory analysis
- Compatible with Seurat and SingleCellExperiment objects

## Installation

```r
# Install from GitHub
devtools::install_github("guokai8/CellHet")

# Load the package
library(CellHet)
```

## Quick Start

Complete workflow with realistic PBMC data (real human gene symbols: CD3D, CD4, CD8A, MS4A1, etc.):

```r
library(CellHet)
library(Seurat)

# Create example PBMC data
pbmc <- createExampleData(
  n_cells = 500,
  n_genes = 200,
  n_cell_types = 4,
  n_conditions = 2,
  seed = 123
)

# Check the data
print(pbmc)
table(pbmc$cell_type, pbmc$condition)

# Run differential expression analysis
deg_results <- findDifferentialGenes(
  object = pbmc,
  group_var = "condition",
  cell_type_var = "cell_type",
  reference_group = "Control",
  logfc_threshold = 0.25,
  p_val_adj_threshold = 0.05,
  workers = 2
)

# View summary
print(deg_results$summary)
head(deg_results$all_degs)

# Create visualizations
heatmap <- visualizeDEGHeatmap(deg_results, show_counts = TRUE)
barplot <- visualizeDEGBarplot(deg_results)
dotplot <- visualizeDEGDotplot(deg_results)
upset <- visualizeDEGUpset(deg_results, by_cell_type = TRUE)

print(heatmap)

# Export results (requires openxlsx)
# exportDEGResults(deg_results, "results.xlsx", split_by = "cell_type")
```

## Core Functions

### 1. Differential Expression Analysis

#### Basic Comparison
```r
# Compare all cell types across conditions
deg_results <- compareDEGs(
  object = pbmc,
  group_var = "condition",
  cell_type_var = "cell_type",
  test_method = "wilcox",
  logfc_threshold = 0.25,
  p_val_adj_threshold = 0.05
)
```

#### Quick Comparison for Specific Cell Types
```r
# Focus on specific cell types
t_cell_results <- quickCompare(
  object = pbmc,
  cell_types = c("CD4_T", "CD8_T"),
  group1 = "Stimulated",
  group2 = "Control",
  return_plot = TRUE
)

print(t_cell_results$results)
print(t_cell_results$plots$CD4_T)  # Volcano plot
```

#### Convert External Results
```r
# Convert Seurat FindMarkers results to CellHet format
deg_list <- list(
  "CD4_T_Treatment_vs_Control" = seurat_deg_results
)

cellhet_deg <- convertToCellHetDEG(
  deg_list = deg_list,
  logfc_threshold = 0.25
)
```

### 2. Visualization

#### Heatmaps
```r
heatmap <- visualizeDEGHeatmap(
  deg_results = deg_results,
  show_counts = TRUE,
  cluster_rows = TRUE,
  up_color = "darkorange",
  down_color = "deepskyblue"
)
```

#### Barplots
```r
barplot <- visualizeDEGBarplot(
  deg_results = deg_results,
  facet_by = "comparison",
  stacked = FALSE
)
```

#### Dotplots
```r
dotplot <- visualizeDEGDotplot(
  deg_results = deg_results,
  use_triangles = TRUE,
  count_display = "inside"
)
```

#### UpSet Plots
```r
upset <- visualizeDEGUpset(
  deg_results = deg_results,
  by_cell_type = TRUE,
  direction = "up",
  highlight_exclusive = TRUE
)
```

#### Multi-Panel Visualization
```r
patterns <- visualizeDEGPatterns(
  deg_results = deg_results,
  plot_types = c("heatmap", "barplot", "dotplot"),
  ncol = 1
)
```

### 3. Analysis Utilities

#### Summarize DEG Results
```r
# Create summary plots
summary_heatmap <- summarizeDEGs(deg_results, plot_type = "heatmap")
summary_barplot <- summarizeDEGs(deg_results, plot_type = "barplot")
summary_dotplot <- summarizeDEGs(deg_results, plot_type = "dotplot")
```

#### Find Shared DEGs
```r
# Analyze gene overlap across comparisons
shared <- findSharedDEGs(
  deg_results = deg_results,
  min_deg_count = 10,
  plot_type = "upset",  # or "venn", "heatmap", "network"
  direction = "both"
)

# Access results
shared$gene_sets      # Gene lists per comparison
shared$overlaps       # Pairwise overlaps
shared$plot           # Visualization
```

#### Compare Cell Type Responses
```r
# Compare responses to a reference
response <- compareResponses(
  deg_results = deg_results,
  reference_cell_type = "CD4_T",
  plot_type = "heatmap",
  similarity_metric = "jaccard"
)

# Visualize responses across conditions
response_viz <- visualizeReferenceComparisons(
  deg_results = deg_results,
  reference_group = "Control",
  plot_type = "all"
)
```

#### Heterogeneity Mapping
```r
het_map <- visualizeHeterogeneity(
  deg_results = deg_results,
  layout_type = "network",
  level = "both"
)
```

### 4. Pathway Enrichment

Requires `richR` package and annotation databases:

```r
# Get annotation data
library(richR)
go_annot <- getGO(organism = "human")

# Run pathway analysis
pathway_results <- runPathwayAnalysis(
  deg_results = deg_results,
  annot_data = go_annot,
  annot_type = "GO",
  ontology = "BP",
  direction = "both"
)

# Visualize pathways
pathway_plot <- visualizePathwaysDotplot(
  pathway_results = pathway_results,
  top_n = 10,
  color_by = "RichFactor",
  flip = TRUE
)

# Alternative: Use richR directly
richr_results <- runRichREnrichment(
  gene_list = deg_genes,
  organism = "human",
  database = "GO"
)
```

### 5. Multi-Modal Integration

```r
# Simulate multi-modal data
multi_data <- simulateMultiModalData(
  n_cells = 800,
  n_modalities = 2,
  n_features = c(2000, 1000),
  n_cell_types = 4,
  n_conditions = 2
)

# Integrate modalities
integration <- integrateMultiModalData(
  data_list = multi_data,
  modality_names = c("RNA", "ATAC"),
  cell_type_var = "cell_type",
  integration_method = "harmony",  # or "canonical_correlation", "mnn", "seurat"
  dim_reduction = "umap",
  n_dims = 30
)

# Visualize integration
viz <- visualizeIntegration(
  integration_results = integration,
  color_by = c("modality", "cell_type"),
  plot_types = c("dimplot", "density")
)

# Evaluate integration quality
eval <- evaluateIntegration(
  integration_results = integration,
  cell_type_var = "cell_type",
  metrics = c("silhouette", "mixing")
)

# Find conserved features
conserved <- findConservedFeatures(
  integration_results = integration,
  modality_names = c("RNA", "ATAC"),
  n_features = 100
)

# Find multi-modal signatures
signatures <- findMultiModalSignatures(
  integration_results = integration,
  cell_type_var = "cell_type",
  n_features_per_modality = 50
)

# Save and load integration models
exportImportIntegrationModel(
  integration_results = integration,
  file_path = "model.rds",
  mode = "export"
)

loaded_model <- exportImportIntegrationModel(
  file_path = "model.rds",
  mode = "import"
)

# Apply to new data
projected <- applyIntegrationModel(
  new_data = new_object,
  integration_model = loaded_model
)
```

### 6. Pseudotime Analysis

Requires `monocle3` package:

```r
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
trajectory <- analyzeDEGTrajectory(
  sce = pbmc,
  trajectory_var = "pseudotime",
  cell_type_var = "cell_type",
  n_bins = 10,
  smooth_method = "loess"
)

# Plot gene trajectories
gene_traj <- plotGeneTrajectory(
  trajectory_results = trajectory,
  genes = c("CD3D", "MS4A1", "NKG7"),
  plot_type = "line"
)

# Compare trajectories between conditions
diff_traj <- diffTrajectoryAnalysis(
  sce = pbmc,
  trajectory_var = "pseudotime",
  group_var = "condition",
  cell_type_var = "cell_type"
)
```

### 7. Data Simulation

```r
# Simulate single-cell data
sim_data <- simulateSingleCellData(
  n_cells = 1000,
  n_genes = 2000,
  n_cell_types = 5,
  n_conditions = 2,
  pct_de_genes = 0.1,
  effect_size = 1.5,
  batch_effect = TRUE,
  add_trajectory = TRUE,
  return_type = "seurat"
)

# Create quick example with real gene names
pbmc_example <- createExampleData(
  n_cells = 500,
  n_genes = 200,
  n_cell_types = 3,
  n_conditions = 2
)
```

### 8. Export and Utilities

```r
# Export DEG results to Excel (requires openxlsx)
exportDEGResults(
  deg_results = deg_results,
  file_path = "results.xlsx",
  split_by = "cell_type",
  include_pathways = TRUE
)

# Export pathway results
exportEnrichmentResults(
  pathway_results = pathway_results,
  file_path = "pathways.xlsx",
  split_by = "cell_type"
)

# Save pathway plot
savePathwayDotplot(
  pathway_plot = pathway_plot,
  file_path = "pathways.pdf",
  width = 10,
  height = 8
)

# Color palettes
colors <- cellhet_colors()              # All colors
colors_10 <- cellhet_colors(10)         # First 10
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

### 9. Interactive Shiny App

```r
# Launch interactive app
runCellHetApp()

# Launch on specific port
runCellHetApp(port = 3838)

# Launch without browser
runCellHetApp(launch.browser = FALSE)
```

## Complete Example Workflow

```r
library(CellHet)
library(Seurat)
library(dplyr)

# 1. Create or load data
pbmc <- createExampleData(n_cells = 500, n_genes = 200, seed = 123)

# 2. Run differential expression
deg_results <- findDifferentialGenes(
  object = pbmc,
  group_var = "condition",
  cell_type_var = "cell_type",
  reference_group = "Control",
  workers = 2
)

# 3. Examine results
print(deg_results$summary)

# View specific cell type
cd4_degs <- deg_results$all_degs %>%
  filter(cell_type == "CD4_T", significant == TRUE) %>%
  arrange(desc(abs(avg_log2FC)))
print(head(cd4_degs, 10))

# 4. Visualize
heatmap <- visualizeDEGHeatmap(deg_results, show_counts = TRUE)
barplot <- visualizeDEGBarplot(deg_results, facet_by = "comparison")
upset <- visualizeDEGUpset(deg_results, by_cell_type = TRUE)

print(heatmap)
print(barplot)
print(upset)

# 5. Analyze shared DEGs
shared <- findSharedDEGs(deg_results, plot_type = "upset", direction = "both")
print(names(shared$gene_sets))
print(head(shared$gene_sets[[1]]))

# 6. Compare cell type responses
response <- compareResponses(deg_results, reference_cell_type = "CD4_T")
response_viz <- visualizeReferenceComparisons(deg_results, reference_group = "Control")

# 7. Export results (requires openxlsx)
# exportDEGResults(deg_results, "pbmc_results.xlsx", split_by = "cell_type")
```

## Function Reference

### Core Analysis Functions
- `compareDEGs()` - Compare DEGs across cell types and conditions
- `findDifferentialGenes()` - Comprehensive DEG analysis pipeline
- `quickCompare()` - Quick comparison for specific cell types
- `convertToCellHetDEG()` - Convert external DEG results

### Visualization Functions
- `visualizeDEGHeatmap()` - DEG count heatmap
- `visualizeDEGBarplot()` - DEG count barplot
- `visualizeDEGDotplot()` - DEG count dotplot
- `visualizeDEGUpset()` - UpSet plot for gene overlaps
- `visualizeDEGPatterns()` - Multi-panel visualization
- `upsetPlot()` - Custom UpSet plot

### Analysis Utilities
- `summarizeDEGs()` - Summary visualizations
- `findSharedDEGs()` - Shared/unique DEG analysis
- `compareResponses()` - Compare cell type responses
- `visualizeReferenceComparisons()` - Reference comparison visualization
- `visualizeHeterogeneity()` - Heterogeneity mapping

### Pathway Analysis
- `runPathwayAnalysis()` - Pathway enrichment analysis
- `runRichREnrichment()` - richR integration
- `visualizePathwaysDotplot()` - Pathway dotplot
- `savePathwayDotplot()` - Save pathway plot
- `exportEnrichmentResults()` - Export pathway results

### Multi-Modal Integration
- `integrateMultiModalData()` - Integrate multiple modalities
- `visualizeIntegration()` - Visualize integration
- `evaluateIntegration()` - Evaluate integration quality
- `optimizeIntegrationParameters()` - Optimize parameters
- `findConservedFeatures()` - Find conserved features
- `findMultiModalSignatures()` - Find multi-modal signatures
- `calculateModalityAgreement()` - Calculate agreement
- `analyzeFeatureImportance()` - Feature importance
- `exportImportIntegrationModel()` - Save/load models
- `applyIntegrationModel()` - Apply to new data
- `projectNewData()` - Project new data

### Pseudotime Analysis
- `calculatePseudotime()` - Calculate pseudotime trajectories
- `plotPseudotime()` - Plot pseudotime
- `getPseudotimeStats()` - Pseudotime statistics
- `analyzeDEGTrajectory()` - Analyze gene trajectories
- `plotGeneTrajectory()` - Plot gene trajectories
- `diffTrajectoryAnalysis()` - Compare trajectories

### Data Simulation
- `simulateSingleCellData()` - Simulate single-cell data
- `simulateMultiModalData()` - Simulate multi-modal data
- `createExampleData()` - Create example PBMC data

### Utilities
- `exportDEGResults()` - Export DEG results to Excel
- `cellhet_colors()` - Get color palette
- `scale_color_cellhet()` - ggplot2 color scale
- `scale_fill_cellhet()` - ggplot2 fill scale
- `runCellHetApp()` - Launch Shiny app

## Citation

If you use CellHet in your research, please cite:

[https://github.com/guokai8/CellHet](https://github.com/guokai8/CellHet)

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

For questions, issues, or feature requests, please open an issue on the [GitHub repository](https://github.com/guokai8/CellHet/issues).
