# CellHet Shiny Application

An interactive web interface for comprehensive single-cell heterogeneity analysis.

## Features

### 1. Data Upload & Overview
- Load example datasets (PBMC)
- Upload your own Seurat or SingleCellExperiment objects (.rds files)
- View data summaries and distribution plots
- Configure metadata variables

### 2. Differential Expression Analysis (DEG)
- Compare gene expression across conditions and cell types
- Set custom parameters (test method, thresholds, min cells)
- Use reference group for automatic comparisons
- Visualizations:
  - Summary barplot of DEG counts
  - Interactive heatmaps
  - Volcano plots
  - UpSet plots for overlap analysis
  - Downloadable data tables

### 3. Pathway Enrichment Analysis
- Perform GO and KEGG enrichment using richR
- Multiple organism support (Human, Mouse)
- Interactive visualizations:
  - Barplots
  - Dotplots
  - Network plots
- Export enrichment results

### 4. Trajectory Analysis
- Analyze gene expression dynamics along pseudotime
- Two analysis types:
  - Basic gene expression trends
  - Differential trajectories between conditions
- Customizable binning and smoothing
- Visualizations:
  - Expression heatmaps
  - Line plots for specific genes
  - Bin summary plots

### 5. Heterogeneity Metrics
- Calculate diversity metrics:
  - Shannon Entropy
  - Simpson Index
  - Gini Coefficient
- Analysis at cell type or group level
- Interactive plots and summary tables

## Launch the App

From R console:

```r
# Load the package
library(CellHet)

# Launch the app
runCellHetApp()

# Or specify custom port
runCellHetApp(port = 3838)
```

## Required Packages

The app requires the following R packages:
- shiny (>= 1.7.0)
- shinythemes
- DT
- ggplot2
- Seurat
- SingleCellExperiment

Optional packages for full functionality:
- richR (for pathway analysis)
- ComplexHeatmap (for advanced heatmaps)
- UpSetR (for UpSet plots)

## Data Input Requirements

### Seurat Objects
- Must be a valid Seurat object (v4+)
- Should contain:
  - Normalized counts (RNA assay)
  - Cell type annotations in metadata
  - Condition/group information in metadata
  - (Optional) Pseudotime values for trajectory analysis

### SingleCellExperiment Objects
- Valid SCE object
- logcounts or counts assay
- colData with cell type and condition information

## Workflow Example

1. **Upload Data**: Load your dataset or use the example PBMC data
2. **Configure Variables**: Select which metadata columns represent cell types and conditions
3. **Run DEG Analysis**:
   - Choose a reference group (e.g., "Control")
   - Set thresholds (log FC, p-value)
   - Run analysis and explore results
4. **Pathway Enrichment**:
   - Select pathway database (GO/KEGG)
   - Run enrichment on DEG results
   - Visualize enriched pathways
5. **Trajectory Analysis** (if pseudotime available):
   - Select pseudotime variable
   - Configure bins and smoothing
   - Visualize gene dynamics
6. **Download Results**: Export plots and tables for publication

## Tips

- Start with example data to familiarize yourself with the interface
- DEG analysis with reference groups simplifies interpretation
- Use at least 10 bins for trajectory analysis for smooth curves
- Download results frequently - data is not saved between sessions
- For large datasets, consider running analyses programmatically first

## Troubleshooting

### App won't launch
- Check that all required packages are installed
- Ensure CellHet package is properly installed
- Try specifying a different port: `runCellHetApp(port = 8080)`

### Analysis fails
- Verify your data has required metadata columns
- Ensure sufficient cells per group (min 3 recommended)
- Check that gene names are valid

### Slow performance
- The app processes data in real-time
- Large datasets (>50,000 cells) may be slow
- Consider filtering to specific cell types first
- Use fewer bins for trajectory analysis

## Contact

For questions, issues, or feature requests:
- GitHub: https://github.com/guokai8/CellHet/issues
- Email: guokai8@gmail.com
