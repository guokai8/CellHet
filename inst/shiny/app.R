#' CellHet Shiny App
#'
#' Interactive interface for CellHet analysis
#'
#' @import shiny
#' @import ggplot2

library(shiny)
library(CellHet)
library(ggplot2)

# Load richR for pathway enrichment analysis
# This ensures all richR dependencies are properly loaded
if (requireNamespace("richR", quietly = TRUE)) {
  suppressPackageStartupMessages(library(richR))
}

# UI Definition
ui <- navbarPage(
  "CellHet: Single-Cell Heterogeneity Analysis",
  theme = shinythemes::shinytheme("flatly"),

  # Tab 1: Data Upload & Overview
  tabPanel("Data Upload",
    sidebarLayout(
      sidebarPanel(width = 3,
        h3("Upload Data"),

        radioButtons("data_source", "Data Source:",
                     choices = c("Example Data" = "example",
                                "Upload RDS" = "upload",
                                "Use Active Seurat Object" = "active"),
                     selected = "example"),

        conditionalPanel(
          condition = "input.data_source == 'upload'",
          fileInput("seurat_file", "Upload Seurat/SCE Object (.rds)",
                   accept = c(".rds", ".RDS"))
        ),

        conditionalPanel(
          condition = "input.data_source == 'example'",
          selectInput("example_data", "Select Example Dataset:",
                     choices = c("PBMC Example (500 cells, 3 cell types)" = "simulated")),
          helpText("Simulated PBMC dataset with realistic gene names (CD4 T, CD8 T, B cells), 3 conditions (Control, Stimulated, Treatment), and pseudotime.")
        ),

        hr(),

        h4("Metadata Variables"),
        uiOutput("group_var_ui"),
        uiOutput("cell_type_var_ui"),

        actionButton("load_data", "Load Data",
                    class = "btn-primary btn-block")
      ),

      mainPanel(
        h3("Data Overview"),
        verbatimTextOutput("data_summary"),

        h4("Cell Type Distribution"),
        plotOutput("cell_type_plot", height = "300px"),

        h4("Group Distribution"),
        plotOutput("group_plot", height = "300px")
      )
    )
  ),

  # Tab 2: DEG Analysis
  tabPanel("DEG Analysis",
    sidebarLayout(
      sidebarPanel(width = 3,
        h3("DEG Parameters"),

        uiOutput("reference_group_ui"),

        numericInput("min_cells_per_group", "Min Cells per Group:",
                    value = 3, min = 1, max = 100),

        selectInput("test_method", "Test Method:",
                   choices = c("Wilcoxon" = "wilcox",
                              "t-test" = "t",
                              "DESeq2" = "deseq2"),
                   selected = "wilcox"),

        numericInput("logfc_threshold", "Log FC Threshold:",
                    value = 0.25, min = 0, max = 5, step = 0.05),

        radioButtons("pval_type", "P-value Type:",
                    choices = c("Adjusted P-value (FDR)" = "adjusted",
                               "Unadjusted P-value" = "unadjusted"),
                    selected = "adjusted"),

        conditionalPanel(
          condition = "input.pval_type == 'adjusted'",
          numericInput("p_val_adj_threshold", "Adjusted P-value Threshold:",
                      value = 0.05, min = 0, max = 1, step = 0.01)
        ),

        conditionalPanel(
          condition = "input.pval_type == 'unadjusted'",
          numericInput("p_val_threshold", "P-value Threshold:",
                      value = 0.05, min = 0, max = 1, step = 0.01)
        ),

        hr(),

        actionButton("run_deg", "Run DEG Analysis",
                    class = "btn-success btn-block"),

        br(),

        conditionalPanel(
          condition = "output.deg_complete",
          selectInput("deg_download_format", "Download Format:",
                     choices = c("Excel (.xlsx)" = "xlsx",
                                "CSV (.csv)" = "csv",
                                "RDS (.rds)" = "rds"),
                     selected = "xlsx"),
          downloadButton("download_deg", "Download Results",
                        class = "btn-info btn-block")
        )
      ),

      mainPanel(
        h3("DEG Results"),

        tabsetPanel(
          tabPanel("Summary",
                  h4("Number of DEGs per Comparison"),
                  plotOutput("deg_summary_plot", height = "400px"),
                  verbatimTextOutput("deg_summary_text")
          ),

          tabPanel("Heatmap",
                  fluidRow(
                    column(8,
                      h4("Top DEGs Heatmap"),
                      uiOutput("heatmap_plot_ui")
                    ),
                    column(4,
                      wellPanel(
                        h5("Plot Settings"),
                        sliderInput("n_top_genes", "Number of Top Genes:",
                                   min = 10, max = 100, value = 50, step = 10),
                        numericInput("heatmap_width", "Figure Width (inches):",
                                    value = 10, min = 4, max = 20, step = 1),
                        numericInput("heatmap_height", "Figure Height (inches):",
                                    value = 8, min = 4, max = 20, step = 1),
                        sliderInput("heatmap_fontsize", "Font Size:",
                                   min = 8, max = 20, value = 12, step = 1),
                        hr(),
                        h5("Color Settings"),
                        colourpicker::colourInput("heatmap_low_color", "Low Value Color:",
                                                 value = "#4DA8A0", showColour = "background"),
                        colourpicker::colourInput("heatmap_mid_color", "Mid Value Color:",
                                                 value = "white", showColour = "background"),
                        colourpicker::colourInput("heatmap_high_color", "High Value Color:",
                                                 value = "#E67E22", showColour = "background"),
                        downloadButton("download_heatmap_pdf", "Download PDF",
                                     class = "btn-info btn-block")
                      )
                    )
                  )
          ),

          tabPanel("Volcano Plot",
                  fluidRow(
                    column(8,
                      h4("Volcano Plot"),
                      uiOutput("volcano_plot_ui")
                    ),
                    column(4,
                      wellPanel(
                        h5("Plot Settings"),
                        uiOutput("comparison_select_ui"),
                        numericInput("volcano_width", "Figure Width (inches):",
                                    value = 8, min = 4, max = 20, step = 1),
                        numericInput("volcano_height", "Figure Height (inches):",
                                    value = 6, min = 4, max = 20, step = 1),
                        sliderInput("volcano_fontsize", "Font Size:",
                                   min = 8, max = 20, value = 12, step = 1),
                        hr(),
                        h5("Color Settings"),
                        colourpicker::colourInput("volcano_up_color", "Up-regulated Color:",
                                                 value = "#E67E22", showColour = "background"),
                        colourpicker::colourInput("volcano_down_color", "Down-regulated Color:",
                                                 value = "#00ACC1", showColour = "background"),
                        colourpicker::colourInput("volcano_ns_color", "Not Significant Color:",
                                                 value = "gray70", showColour = "background"),
                        downloadButton("download_volcano_pdf", "Download PDF",
                                     class = "btn-info btn-block")
                      )
                    )
                  )
          ),

          tabPanel("UpSet Plot",
                  fluidRow(
                    column(8,
                      h4("DEG Overlap (UpSet Plot)"),
                      uiOutput("upset_plot_ui")
                    ),
                    column(4,
                      wellPanel(
                        h5("Plot Settings"),
                        numericInput("upset_width", "Figure Width (inches):",
                                    value = 10, min = 4, max = 20, step = 1),
                        numericInput("upset_height", "Figure Height (inches):",
                                    value = 6, min = 4, max = 20, step = 1),
                        sliderInput("upset_fontsize", "Font Size:",
                                   min = 8, max = 20, value = 12, step = 1),
                        hr(),
                        h5("Color Settings"),
                        colourpicker::colourInput("upset_bar_color", "Intersection Bar Color:",
                                                 value = "#2670B8", showColour = "background"),
                        colourpicker::colourInput("upset_point_color", "Point/Line Color:",
                                                 value = "black", showColour = "background"),
                        checkboxInput("upset_use_palette", "Use Color Palette for Sets", value = TRUE),
                        downloadButton("download_upset_pdf", "Download PDF",
                                     class = "btn-info btn-block")
                      )
                    )
                  )
          ),

          tabPanel("Data Table",
                  h4("DEG Table"),
                  DT::dataTableOutput("deg_table")
          )
        )
      )
    )
  ),

  # Tab 3: Pathway Analysis
  tabPanel("Pathway Analysis",
    sidebarLayout(
      sidebarPanel(width = 3,
        h3("Pathway Parameters"),

        selectInput("pathway_type", "Pathway Database:",
                   choices = c("GO Biological Process" = "GO_BP",
                              "GO Molecular Function" = "GO_MF",
                              "GO Cellular Component" = "GO_CC",
                              "KEGG" = "KEGG"),
                   selected = "GO_BP"),

        selectInput("pathway_organism", "Organism:",
                   choices = c("Human" = "human",
                              "Mouse" = "mouse"),
                   selected = "human"),

        numericInput("pathway_pval", "P-value Threshold:",
                    value = 0.05, min = 0, max = 1, step = 0.01),

        numericInput("top_pathways", "Top Pathways to Show:",
                    value = 20, min = 5, max = 50, step = 5),

        hr(),

        actionButton("run_pathway", "Run Pathway Analysis",
                    class = "btn-success btn-block"),

        br(),

        conditionalPanel(
          condition = "output.pathway_complete",
          downloadButton("download_pathway", "Download Results",
                        class = "btn-info btn-block")
        )
      ),

      mainPanel(
        h3("Pathway Enrichment Results"),

        tabsetPanel(
          tabPanel("Barplot",
                  fluidRow(
                    column(8,
                      h4("Pathway Barplot"),
                      uiOutput("pathway_barplot_ui")
                    ),
                    column(4,
                      wellPanel(
                        h5("Plot Settings"),
                        uiOutput("pathway_comparison_ui"),
                        numericInput("pathway_bar_width", "Figure Width (inches):",
                                    value = 10, min = 4, max = 20, step = 1),
                        numericInput("pathway_bar_height", "Figure Height (inches):",
                                    value = 8, min = 4, max = 20, step = 1),
                        sliderInput("pathway_bar_fontsize", "Font Size:",
                                   min = 8, max = 20, value = 12, step = 1),
                        hr(),
                        h5("Color Settings"),
                        colourpicker::colourInput("pathway_bar_fill_color", "Bar Fill Color:",
                                                 value = "#2670B8", showColour = "background"),
                        downloadButton("download_pathway_bar_pdf", "Download PDF",
                                     class = "btn-info btn-block")
                      )
                    )
                  )
          ),

          tabPanel("Dotplot",
                  fluidRow(
                    column(8,
                      h4("Pathway Dotplot"),
                      uiOutput("pathway_dotplot_ui")
                    ),
                    column(4,
                      wellPanel(
                        h5("Plot Settings"),
                        numericInput("pathway_dot_width", "Figure Width (inches):",
                                    value = 10, min = 4, max = 20, step = 1),
                        numericInput("pathway_dot_height", "Figure Height (inches):",
                                    value = 8, min = 4, max = 20, step = 1),
                        sliderInput("pathway_dot_fontsize", "Font Size:",
                                   min = 8, max = 20, value = 12, step = 1),
                        hr(),
                        h5("Color Settings"),
                        colourpicker::colourInput("pathway_dot_low_color", "Low P-value Color:",
                                                 value = "#E67E22", showColour = "background"),
                        colourpicker::colourInput("pathway_dot_high_color", "High P-value Color:",
                                                 value = "#2670B8", showColour = "background"),
                        downloadButton("download_pathway_dot_pdf", "Download PDF",
                                     class = "btn-info btn-block")
                      )
                    )
                  )
          ),

          tabPanel("Network",
                  fluidRow(
                    column(8,
                      h4("Pathway Network"),
                      uiOutput("pathway_network_ui")
                    ),
                    column(4,
                      wellPanel(
                        h5("Plot Settings"),
                        numericInput("pathway_net_width", "Figure Width (inches):",
                                    value = 12, min = 4, max = 20, step = 1),
                        numericInput("pathway_net_height", "Figure Height (inches):",
                                    value = 10, min = 4, max = 20, step = 1),
                        sliderInput("pathway_net_fontsize", "Font Size:",
                                   min = 8, max = 20, value = 12, step = 1),
                        hr(),
                        helpText("Note: Network plot colors are determined by richR and cannot be customized."),
                        downloadButton("download_pathway_net_pdf", "Download PDF",
                                     class = "btn-info btn-block")
                      )
                    )
                  )
          ),

          tabPanel("Table",
                  DT::dataTableOutput("pathway_table")
          )
        )
      )
    )
  ),

  # Tab 4: Trajectory Analysis
  tabPanel("Trajectory Analysis",
    sidebarLayout(
      sidebarPanel(width = 3,
        h3("Trajectory Parameters"),

        uiOutput("trajectory_var_ui"),

        selectInput("trajectory_type", "Analysis Type:",
                   choices = c("Gene Expression Trends" = "basic",
                              "Differential Trajectories" = "differential"),
                   selected = "basic"),

        numericInput("n_bins", "Number of Bins:",
                    value = 10, min = 3, max = 20, step = 1),

        numericInput("min_cells_per_bin", "Min Cells per Bin:",
                    value = 5, min = 1, max = 50, step = 1),

        selectInput("smooth_method", "Smoothing Method:",
                   choices = c("LOESS" = "loess",
                              "GAM" = "gam",
                              "None" = "none"),
                   selected = "loess"),

        numericInput("n_top_trajectory_genes", "Top Variable Genes:",
                    value = 50, min = 10, max = 200, step = 10),

        hr(),

        actionButton("run_trajectory", "Run Trajectory Analysis",
                    class = "btn-success btn-block")
      ),

      mainPanel(
        h3("Trajectory Analysis Results"),

        tabsetPanel(
          tabPanel("Heatmap",
                  fluidRow(
                    column(9,
                      plotOutput("trajectory_heatmap", height = "600px")
                    ),
                    column(3,
                      wellPanel(
                        h5("Color Settings"),
                        colourpicker::colourInput("traj_low_color", "Low Expression Color:",
                                                 value = "deepskyblue", showColour = "background"),
                        colourpicker::colourInput("traj_mid_color", "Mid Expression Color:",
                                                 value = "white", showColour = "background"),
                        colourpicker::colourInput("traj_high_color", "High Expression Color:",
                                                 value = "darkorange", showColour = "background")
                      )
                    )
                  )
          ),

          tabPanel("Line Plots",
                  uiOutput("gene_select_ui"),
                  plotOutput("trajectory_lineplot", height = "500px")
          ),

          tabPanel("Summary",
                  verbatimTextOutput("trajectory_summary"),
                  plotOutput("trajectory_bins_plot", height = "300px")
          )
        )
      )
    )
  ),

  # Tab 5: Heterogeneity Analysis
  tabPanel("Heterogeneity",
    sidebarLayout(
      sidebarPanel(width = 3,
        h3("Heterogeneity Analysis"),

        helpText("Visualize heterogeneity across cell types based on DEG similarity"),

        selectInput("het_layout", "Visualization Type:",
                   choices = c("Heatmap" = "heatmap",
                              "Network" = "network",
                              "UMAP" = "umap",
                              "Circular" = "circular"),
                   selected = "heatmap"),

        selectInput("het_level", "Analysis Level:",
                   choices = c("Gene" = "gene",
                              "Pathway" = "pathway",
                              "Both" = "both"),
                   selected = "gene"),

        selectInput("het_metric", "Similarity Metric:",
                   choices = c("Jaccard" = "jaccard",
                              "Overlap" = "overlap",
                              "Correlation" = "correlation"),
                   selected = "jaccard"),

        hr(),

        actionButton("run_het", "Calculate Heterogeneity",
                    class = "btn-success btn-block")
      ),

      mainPanel(
        h3("Heterogeneity Results"),

        tabsetPanel(
          tabPanel("Plot",
                  plotOutput("het_plot", height = "600px")
          ),

          tabPanel("Matrix",
                  h4("Heterogeneity Matrix"),
                  DT::dataTableOutput("het_matrix_table"),
                  br(),
                  h4("Similarity Matrix"),
                  DT::dataTableOutput("het_similarity_table")
          )
        )
      )
    )
  ),

  # Tab 6: Help & About
  tabPanel("Help",
    fluidRow(
      column(12,
        h2("CellHet: Single-Cell Heterogeneity Analysis"),

        h3("Overview"),
        p("CellHet is a comprehensive R package for analyzing cellular heterogeneity in single-cell RNA-seq data."),

        h3("Key Features"),
        tags$ul(
          tags$li(strong("DEG Analysis:"), " Identify differentially expressed genes across conditions and cell types"),
          tags$li(strong("Pathway Enrichment:"), " Perform GO/KEGG enrichment analysis using richR"),
          tags$li(strong("Trajectory Analysis:"), " Analyze gene expression dynamics along pseudotime"),
          tags$li(strong("Heterogeneity Metrics:"), " Quantify cellular diversity and heterogeneity")
        ),

        h3("Workflow"),
        tags$ol(
          tags$li(strong("Upload Data:"), " Load your Seurat or SingleCellExperiment object"),
          tags$li(strong("Configure Parameters:"), " Set analysis parameters for each module"),
          tags$li(strong("Run Analysis:"), " Execute analyses and explore interactive visualizations"),
          tags$li(strong("Download Results:"), " Export results and plots for further use")
        ),

        h3("Citation"),
        p("If you use CellHet in your research, please cite:"),
        tags$pre("CellHet: A Comprehensive Package for Single-Cell Heterogeneity Analysis"),

        h3("Support"),
        p("For questions, issues, or feature requests, please visit our GitHub repository."),

        hr(),

        h4("Package Information"),
        verbatimTextOutput("session_info")
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {

  # Reactive values to store data
  values <- reactiveValues(
    seurat_obj = NULL,
    metadata = NULL,
    deg_results = NULL,
    pathway_results = NULL,
    trajectory_results = NULL,
    het_results = NULL
  )

  # ==================== Data Upload ====================

  # Render group variable selector
  output$group_var_ui <- renderUI({
    req(values$metadata)
    selectInput("group_var", "Group Variable:",
               choices = colnames(values$metadata),
               selected = "condition")
  })

  # Render cell type variable selector
  output$cell_type_var_ui <- renderUI({
    req(values$metadata)
    selectInput("cell_type_var", "Cell Type Variable:",
               choices = colnames(values$metadata),
               selected = "cell_type")
  })

  # Load data
  observeEvent(input$load_data, {

    withProgress(message = 'Loading data...', value = 0, {

      if (input$data_source == "example") {
        # Load example data
        incProgress(0.3, detail = "Creating example dataset")

        tryCatch({
          values$seurat_obj <- CellHet::createExampleData(
            n_cells = 500,
            n_genes = 200,
            n_cell_types = 3,
            n_conditions = 3,
            seed = 123
          )
          incProgress(0.5, detail = "Example data created")

        }, error = function(e) {
          showNotification(
            paste("Error creating example data:", e$message),
            type = "error",
            duration = 10
          )
          return(NULL)
        })

      } else if (input$data_source == "upload") {
        # Load uploaded file
        req(input$seurat_file)
        incProgress(0.3, detail = "Reading uploaded file")
        values$seurat_obj <- readRDS(input$seurat_file$datapath)
      }

      incProgress(0.6, detail = "Extracting metadata")

      # Extract metadata
      if (inherits(values$seurat_obj, "Seurat")) {
        values$metadata <- values$seurat_obj@meta.data
      } else if (inherits(values$seurat_obj, "SingleCellExperiment")) {
        values$metadata <- as.data.frame(SingleCellExperiment::colData(values$seurat_obj))
      }

      incProgress(1, detail = "Done!")
    })

    showNotification("Data loaded successfully!", type = "message")
  })

  # Data summary
  output$data_summary <- renderPrint({
    req(values$seurat_obj)

    cat("Object Type:", class(values$seurat_obj)[1], "\n")
    cat("Number of Cells:", ncol(values$seurat_obj), "\n")
    cat("Number of Features:", nrow(values$seurat_obj), "\n\n")

    cat("Available Metadata Variables:\n")
    print(colnames(values$metadata))
  })

  # Cell type distribution plot
  output$cell_type_plot <- renderPlot({
    req(values$metadata, input$cell_type_var)

    df <- data.frame(table(values$metadata[[input$cell_type_var]]))
    colnames(df) <- c("CellType", "Count")

    ggplot(df, aes(x = reorder(CellType, -Count), y = Count, fill = CellType)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = CellHet::cellhet_colors()) +
      theme_minimal() +
      labs(x = "Cell Type", y = "Number of Cells") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
  })

  # Group distribution plot
  output$group_plot <- renderPlot({
    req(values$metadata, input$group_var)

    df <- data.frame(table(values$metadata[[input$group_var]]))
    colnames(df) <- c("Group", "Count")

    ggplot(df, aes(x = Group, y = Count, fill = Group)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = CellHet::cellhet_colors()) +
      theme_minimal() +
      labs(x = "Group", y = "Number of Cells") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
  })

  # ==================== DEG Analysis ====================

  # Reference group selector
  output$reference_group_ui <- renderUI({
    req(values$metadata, input$group_var)

    groups <- unique(values$metadata[[input$group_var]])
    selectInput("reference_group", "Reference Group (optional):",
               choices = c("None" = "", as.character(groups)),
               selected = "")
  })

  # Run DEG analysis
  observeEvent(input$run_deg, {
    req(values$seurat_obj, input$group_var, input$cell_type_var)

    withProgress(message = 'Running DEG analysis...', value = 0, {

      incProgress(0.2, detail = "Preparing data")

      ref_group <- if (input$reference_group == "") NULL else input$reference_group

      incProgress(0.4, detail = "Running comparisons")

      tryCatch({
        # Determine which p-value threshold to use
        p_val_threshold_param <- NULL
        p_val_adj_threshold_param <- 0.05

        if (input$pval_type == "unadjusted") {
          p_val_threshold_param <- input$p_val_threshold
        } else {
          p_val_adj_threshold_param <- input$p_val_adj_threshold
        }

        values$deg_results <- compareDEGs(
          object = values$seurat_obj,
          group_var = input$group_var,
          cell_type_var = input$cell_type_var,
          reference_group = ref_group,
          min_cells_per_group = input$min_cells_per_group,
          test_method = input$test_method,
          logfc_threshold = input$logfc_threshold,
          p_val_threshold = p_val_threshold_param,
          p_val_adj_threshold = p_val_adj_threshold_param
        )

        incProgress(1, detail = "Done!")
        showNotification("DEG analysis completed!", type = "message")

      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    })
  })

  # DEG complete flag
  output$deg_complete <- reactive({
    !is.null(values$deg_results)
  })
  outputOptions(output, "deg_complete", suspendWhenHidden = FALSE)

  # DEG summary plot
  output$deg_summary_plot <- renderPlot({
    req(values$deg_results)

    # Use summary data from compareDEGs results
    if (!is.null(values$deg_results$summary) && nrow(values$deg_results$summary) > 0) {
      # Use the summary table
      plot_data <- values$deg_results$summary
      plot_data$comparison_label <- paste(plot_data$cell_type, plot_data$comparison, sep = " - ")

      ggplot(plot_data, aes(x = reorder(comparison_label, -n_significant),
                           y = n_significant,
                           fill = cell_type)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = CellHet::cellhet_colors()) +
        theme_minimal() +
        labs(x = "Comparison", y = "Number of Significant DEGs", fill = "Cell Type") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    } else {
      # Fallback: empty plot with message
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                label = "No DEG results to display",
                size = 6) +
        theme_void()
    }
  })

  # DEG summary text
  output$deg_summary_text <- renderPrint({
    req(values$deg_results)

    cat("Total Comparisons:", length(values$deg_results$degs), "\n")
    cat("Total DEGs:", sum(values$deg_results$summary$n_significant, na.rm = TRUE), "\n")
    cat("Up-regulated:", sum(values$deg_results$summary$n_up, na.rm = TRUE), "\n")
    cat("Down-regulated:", sum(values$deg_results$summary$n_down, na.rm = TRUE), "\n")
  })

  # DEG heatmap - reactive UI for dynamic height and width
  output$heatmap_plot_ui <- renderUI({
    width_px <- if (!is.null(input$heatmap_width)) input$heatmap_width * 100 else 1000
    height_px <- if (!is.null(input$heatmap_height)) input$heatmap_height * 100 else 800
    plotOutput("deg_heatmap", width = paste0(width_px, "px"), height = paste0(height_px, "px"))
  })

  # DEG heatmap
  output$deg_heatmap <- renderPlot({
    req(values$deg_results)

    tryCatch({
      # Get settings from input
      n_genes <- input$n_top_genes
      font_size <- if (!is.null(input$heatmap_fontsize)) input$heatmap_fontsize else 12

      # Get color settings
      low_col <- if (!is.null(input$heatmap_low_color)) input$heatmap_low_color else "#4DA8A0"
      mid_col <- if (!is.null(input$heatmap_mid_color)) input$heatmap_mid_color else "white"
      high_col <- if (!is.null(input$heatmap_high_color)) input$heatmap_high_color else "#E67E22"

      # Create heatmap using visualizeDEGHeatmap with font size and colors
      CellHet::visualizeDEGHeatmap(
        deg_results = values$deg_results,
        comparison_col = "comparison",
        cell_type_col = "cell_type",
        show_counts = TRUE,
        cluster_rows = TRUE,
        low_color = low_col,
        mid_color = mid_col,
        high_color = high_col,
        title = paste("Top", n_genes, "DEGs Across Cell Types and Comparisons"),
        title_size = font_size + 2,
        count_size = font_size / 4,
        top_label_size = font_size / 3
      )
    }, error = function(e) {
      # If heatmap fails, show empty plot with error message
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                label = paste("Error creating heatmap:", e$message),
                size = 5) +
        theme_void()
    })
  })

  # Comparison selector for volcano plot
  output$comparison_select_ui <- renderUI({
    req(values$deg_results)

    # Get all available comparisons
    comparison_keys <- names(values$deg_results$degs)

    selectInput("selected_comparison", "Select Comparison:",
               choices = comparison_keys,
               selected = comparison_keys[1])
  })

  # Volcano plot - reactive UI for dynamic height and width
  output$volcano_plot_ui <- renderUI({
    width_px <- if (!is.null(input$volcano_width)) input$volcano_width * 100 else 800
    height_px <- if (!is.null(input$volcano_height)) input$volcano_height * 100 else 600
    plotOutput("volcano_plot", width = paste0(width_px, "px"), height = paste0(height_px, "px"))
  })

  # Volcano plot
  output$volcano_plot <- renderPlot({
    req(values$deg_results, input$selected_comparison)

    tryCatch({
      # Get the selected comparison data
      deg_data <- values$deg_results$degs[[input$selected_comparison]]

      if (is.null(deg_data) || nrow(deg_data) == 0) {
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                  label = "No data available for selected comparison",
                  size = 6) +
          theme_void()
      } else {
        # Create volcano plot data
        deg_data$neg_log10_p <- -log10(deg_data$p_val_adj)
        deg_data$color <- "Not Significant"
        deg_data$color[deg_data$significant & deg_data$avg_log2FC > 0] <- "Up-regulated"
        deg_data$color[deg_data$significant & deg_data$avg_log2FC < 0] <- "Down-regulated"

        # Label top genes
        deg_data$label <- ""
        if (sum(deg_data$significant) > 0) {
          sig_data <- deg_data[deg_data$significant, ]
          sig_data <- sig_data[order(abs(sig_data$avg_log2FC), decreasing = TRUE), ]
          top_genes <- head(sig_data$gene, 10)
          deg_data$label[deg_data$gene %in% top_genes] <- deg_data$gene[deg_data$gene %in% top_genes]
        }

        # Determine threshold to display
        p_threshold <- if (!is.null(input$p_val_threshold) && input$pval_type == "unadjusted") {
          -log10(input$p_val_threshold)
        } else {
          -log10(input$p_val_adj_threshold)
        }

        # Get font size setting
        font_size <- if (!is.null(input$volcano_fontsize)) input$volcano_fontsize else 12

        # Get color settings
        up_color <- if (!is.null(input$volcano_up_color)) input$volcano_up_color else "#E67E22"
        down_color <- if (!is.null(input$volcano_down_color)) input$volcano_down_color else "#00ACC1"
        ns_color <- if (!is.null(input$volcano_ns_color)) input$volcano_ns_color else "gray70"

        # Create the plot
        p <- ggplot(deg_data, aes(x = avg_log2FC, y = neg_log10_p, color = color)) +
          geom_point(alpha = 0.6, size = 2) +
          scale_color_manual(values = c("Up-regulated" = up_color,
                                       "Down-regulated" = down_color,
                                       "Not Significant" = ns_color)) +
          geom_hline(yintercept = p_threshold, linetype = "dashed", color = "red") +
          geom_vline(xintercept = c(-input$logfc_threshold, input$logfc_threshold),
                    linetype = "dashed", color = "red") +
          theme_minimal(base_size = font_size) +
          labs(title = paste("Volcano Plot:", input$selected_comparison),
               x = "Log2 Fold Change",
               y = "-Log10 Adjusted P-value",
               color = "Regulation")

        # Add labels if there are any
        if (any(deg_data$label != "")) {
          p <- p + ggrepel::geom_text_repel(
            aes(label = label),
            size = font_size / 4,
            max.overlaps = 20
          )
        }

        p
      }
    }, error = function(e) {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                label = paste("Error creating volcano plot:", e$message),
                size = 5) +
        theme_void()
    })
  })

  # UpSet plot - reactive UI for dynamic height and width
  output$upset_plot_ui <- renderUI({
    width_px <- if (!is.null(input$upset_width)) input$upset_width * 100 else 1000
    height_px <- if (!is.null(input$upset_height)) input$upset_height * 100 else 600
    plotOutput("upset_plot", width = paste0(width_px, "px"), height = paste0(height_px, "px"))
  })

  # UpSet plot
  output$upset_plot <- renderPlot({
    req(values$deg_results)

    tryCatch({
      # Get font size setting
      font_size <- if (!is.null(input$upset_fontsize)) input$upset_fontsize else 12

      # Get color settings
      bar_color <- if (!is.null(input$upset_bar_color)) input$upset_bar_color else "#2670B8"
      point_color <- if (!is.null(input$upset_point_color)) input$upset_point_color else "black"
      use_palette <- if (!is.null(input$upset_use_palette)) input$upset_use_palette else TRUE

      # Set sets_bar_colors based on checkbox
      sets_colors <- if (use_palette) NULL else bar_color

      # Create UpSet plot using visualizeDEGUpset with font size and colors
      CellHet::visualizeDEGUpset(
        deg_results = values$deg_results,
        by_cell_type = TRUE,
        direction = "both",
        min_size = 3,
        bar_width = 0.6,
        main_bar_color = bar_color,
        intersection_color = point_color,
        point_outline_color = point_color,
        sets_bar_colors = sets_colors,
        title = "DEG Overlap Across Cell Types",
        text_size = font_size,
        set_text_size = font_size,
        set_label_size = font_size / 4,
        intersection_label_size = font_size / 4
      )
    }, error = function(e) {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                label = paste("Error creating UpSet plot:", e$message),
                size = 5) +
        theme_void()
    })
  })

  # DEG data table
  output$deg_table <- DT::renderDataTable({
    req(values$deg_results)

    # Combine all DEG results into one table
    all_degs <- data.frame()

    for (result_key in names(values$deg_results$degs)) {
      # Extract cell type and comparison from key
      key_parts <- strsplit(result_key, "__")[[1]]
      cell_type <- key_parts[1]
      comparison <- key_parts[2]

      # Get DEG results
      degs <- values$deg_results$degs[[result_key]]

      # Filter for significant DEGs that pass both thresholds
      if ("passes_lfc_threshold" %in% colnames(degs)) {
        degs_filtered <- degs[degs$significant & degs$passes_lfc_threshold, ]
      } else {
        degs_filtered <- degs[degs$significant, ]
      }

      if (nrow(degs_filtered) > 0) {
        # Add metadata columns
        degs_filtered$cell_type <- cell_type
        degs_filtered$comparison <- comparison

        # Append to combined table
        all_degs <- rbind(all_degs, degs_filtered)
      }
    }

    if (nrow(all_degs) > 0) {
      # Select and reorder columns for display
      display_cols <- c("gene", "cell_type", "comparison", "avg_log2FC", "p_val", "p_val_adj", "direction")
      display_cols <- display_cols[display_cols %in% colnames(all_degs)]

      all_degs <- all_degs[, display_cols]

      # Round numeric columns
      all_degs$avg_log2FC <- round(all_degs$avg_log2FC, 3)
      all_degs$p_val <- formatC(all_degs$p_val, format = "e", digits = 2)
      all_degs$p_val_adj <- formatC(all_degs$p_val_adj, format = "e", digits = 2)

      DT::datatable(all_degs,
                    options = list(pageLength = 25, scrollX = TRUE),
                    rownames = FALSE,
                    filter = "top")
    } else {
      DT::datatable(data.frame(Message = "No significant DEGs found"))
    }
  })

  # Download DEG results handler
  output$download_deg <- downloadHandler(
    filename = function() {
      format_ext <- switch(input$deg_download_format,
                          "xlsx" = ".xlsx",
                          "csv" = ".csv",
                          "rds" = ".rds",
                          ".xlsx")
      paste0("CellHet_DEG_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), format_ext)
    },
    content = function(file) {
      req(values$deg_results, input$deg_download_format)

      withProgress(message = 'Preparing download...', value = 0, {

        incProgress(0.3, detail = "Formatting results")

        tryCatch({
          if (input$deg_download_format == "xlsx") {
            # Export to Excel with multiple sheets
            CellHet::exportDEGResults(
              deg_results = values$deg_results,
              file_path = file,
              split_by = "cell_type",
              include_pathways = FALSE
            )

          } else if (input$deg_download_format == "csv") {
            # Combine all DEG results into one CSV
            all_degs <- data.frame()

            for (result_key in names(values$deg_results$degs)) {
              key_parts <- strsplit(result_key, "__")[[1]]
              cell_type <- key_parts[1]
              comparison <- key_parts[2]

              degs <- values$deg_results$degs[[result_key]]
              degs$cell_type <- cell_type
              degs$comparison <- comparison

              all_degs <- rbind(all_degs, degs)
            }

            write.csv(all_degs, file, row.names = FALSE)

          } else if (input$deg_download_format == "rds") {
            # Save entire results object as RDS
            saveRDS(values$deg_results, file)
          }

          incProgress(1, detail = "Done!")

        }, error = function(e) {
          showNotification(
            paste("Error creating download:", e$message),
            type = "error",
            duration = 10
          )
        })
      })
    }
  )

  # ==================== Pathway Analysis ====================

  # Run pathway analysis
  observeEvent(input$run_pathway, {
    req(values$deg_results)

    withProgress(message = 'Running pathway enrichment...', value = 0, {

      incProgress(0.2, detail = "Preparing gene lists")

      tryCatch({
        # Determine analysis type from pathway_type
        analysis_type <- if (grepl("GO", input$pathway_type)) {
          "GO"
        } else if (input$pathway_type == "KEGG") {
          "KEGG"
        } else {
          "GO"
        }

        # Determine ontology for GO
        ontology <- if (input$pathway_type == "GO_BP") {
          "BP"
        } else if (input$pathway_type == "GO_MF") {
          "MF"
        } else if (input$pathway_type == "GO_CC") {
          "CC"
        } else {
          "BP"
        }

        incProgress(0.4, detail = "Running enrichment analysis")

        # Run enrichment using runRichREnrichment
        values$pathway_results <- CellHet::runRichREnrichment(
          deg_results = values$deg_results,
          analysis_type = analysis_type,
          organism = input$pathway_organism,
          ontology = ontology,
          pvalue = input$pathway_pval,
          qvalue = input$pathway_pval,
          use_significant_only = TRUE,
          separate_direction = FALSE,
          visualize = TRUE,
          plot_type = "all",
          top_terms = input$top_pathways,
          return_plots = TRUE
        )

        incProgress(1, detail = "Done!")
        showNotification("Pathway enrichment completed!", type = "message")

      }, error = function(e) {
        showNotification(paste("Error in pathway analysis:", e$message), type = "error", duration = 10)
      })
    })
  })

  # Pathway complete flag
  output$pathway_complete <- reactive({
    !is.null(values$pathway_results)
  })
  outputOptions(output, "pathway_complete", suspendWhenHidden = FALSE)

  # Pathway comparison selector
  output$pathway_comparison_ui <- renderUI({
    req(values$pathway_results)

    if (!is.null(values$pathway_results$enrichment_results)) {
      # Get available cell types/comparisons
      if ("cell_type" %in% names(values$pathway_results$enrichment_results)) {
        comparison_keys <- names(values$pathway_results$enrichment_results)
        selectInput("selected_pathway_comparison", "Select Cell Type/Comparison:",
                   choices = comparison_keys,
                   selected = comparison_keys[1])
      }
    }
  })

  # Pathway barplot - reactive UI for dynamic height and width
  output$pathway_barplot_ui <- renderUI({
    width_px <- if (!is.null(input$pathway_bar_width)) input$pathway_bar_width * 100 else 1000
    height_px <- if (!is.null(input$pathway_bar_height)) input$pathway_bar_height * 100 else 800
    plotOutput("pathway_barplot", width = paste0(width_px, "px"), height = paste0(height_px, "px"))
  })

  # Pathway barplot
  output$pathway_barplot <- renderPlot({
    req(values$pathway_results)

    tryCatch({
      # Check if plots are available
      if (!is.null(values$pathway_results$plots) && "bar" %in% names(values$pathway_results$plots)) {
        # Apply font size and color to richR plot
        font_size <- if (!is.null(input$pathway_bar_fontsize)) input$pathway_bar_fontsize else 12
        fill_color <- if (!is.null(input$pathway_bar_fill_color)) input$pathway_bar_fill_color else "#2670B8"

        # Modify the plot with custom colors
        values$pathway_results$plots$bar +
          scale_fill_gradient(low = "white", high = fill_color) +
          theme(text = element_text(size = font_size),
                axis.text = element_text(size = font_size),
                axis.title = element_text(size = font_size),
                plot.title = element_text(size = font_size + 2))
      } else {
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                  label = "No enrichment results found. Try adjusting p-value threshold.",
                  size = 6) +
          theme_void()
      }
    }, error = function(e) {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                label = paste("Error creating barplot:", e$message),
                size = 5) +
        theme_void()
    })
  })

  # Pathway dotplot - reactive UI for dynamic height and width
  output$pathway_dotplot_ui <- renderUI({
    width_px <- if (!is.null(input$pathway_dot_width)) input$pathway_dot_width * 100 else 1000
    height_px <- if (!is.null(input$pathway_dot_height)) input$pathway_dot_height * 100 else 800
    plotOutput("pathway_dotplot", width = paste0(width_px, "px"), height = paste0(height_px, "px"))
  })

  # Pathway dotplot
  output$pathway_dotplot <- renderPlot({
    req(values$pathway_results)

    tryCatch({
      if (!is.null(values$pathway_results$plots) && "dot" %in% names(values$pathway_results$plots)) {
        # Apply font size and color to richR plot
        font_size <- if (!is.null(input$pathway_dot_fontsize)) input$pathway_dot_fontsize else 12
        low_color <- if (!is.null(input$pathway_dot_low_color)) input$pathway_dot_low_color else "#E67E22"
        high_color <- if (!is.null(input$pathway_dot_high_color)) input$pathway_dot_high_color else "#2670B8"

        # Modify the plot with custom colors
        values$pathway_results$plots$dot +
          scale_color_gradient(low = low_color, high = high_color) +
          theme(text = element_text(size = font_size),
                axis.text = element_text(size = font_size),
                axis.title = element_text(size = font_size),
                plot.title = element_text(size = font_size + 2),
                legend.text = element_text(size = font_size))
      } else {
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                  label = "No enrichment results found. Try adjusting p-value threshold.",
                  size = 6) +
          theme_void()
      }
    }, error = function(e) {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                label = paste("Error creating dotplot:", e$message),
                size = 5) +
        theme_void()
    })
  })

  # Pathway network plot - reactive UI for dynamic height and width
  output$pathway_network_ui <- renderUI({
    width_px <- if (!is.null(input$pathway_net_width)) input$pathway_net_width * 100 else 1200
    height_px <- if (!is.null(input$pathway_net_height)) input$pathway_net_height * 100 else 1000
    plotOutput("pathway_network", width = paste0(width_px, "px"), height = paste0(height_px, "px"))
  })

  # Pathway network plot
  output$pathway_network <- renderPlot({
    req(values$pathway_results)

    tryCatch({
      if (!is.null(values$pathway_results$plots) && "network" %in% names(values$pathway_results$plots)) {
        # Display the richR network plot as-is
        # Note: richR network plots may not support ggplot2 layer modifications
        print(values$pathway_results$plots$network)
      } else {
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                  label = "Network plot not available. Requires sufficient enriched pathways.",
                  size = 6) +
          theme_void()
      }
    }, error = function(e) {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                label = paste("Error creating network plot:", e$message),
                size = 5) +
        theme_void()
    })
  })

  # Pathway data table
  output$pathway_table <- DT::renderDataTable({
    req(values$pathway_results)

    if (!is.null(values$pathway_results$enrichment_results)) {
      # Combine all enrichment results
      all_enrichment <- data.frame()

      # Check the structure of enrichment_results
      if (is.list(values$pathway_results$enrichment_results)) {
        for (key in names(values$pathway_results$enrichment_results)) {
          result_df <- values$pathway_results$enrichment_results[[key]]
          if (is.data.frame(result_df) && nrow(result_df) > 0) {
            result_df$source <- key
            all_enrichment <- rbind(all_enrichment, result_df)
          }
        }
      } else if (is.data.frame(values$pathway_results$enrichment_results)) {
        all_enrichment <- values$pathway_results$enrichment_results
      }

      if (nrow(all_enrichment) > 0) {
        # Select ALL richR columns for display
        # Core columns (common to GO and KEGG)
        display_cols <- c("Annot", "Term", "Annotated", "Significant", "RichFactor",
                         "FoldEnrichment", "zscore", "p.value", "p.adjust", "GeneID")

        # Add KEGG-specific hierarchical columns if present
        kegg_cols <- c("ko", "Level3", "Level2", "Level1")
        for (col in kegg_cols) {
          if (col %in% colnames(all_enrichment)) {
            display_cols <- c(display_cols, col)
          }
        }

        # Add source column if it exists
        if ("source" %in% colnames(all_enrichment)) {
          display_cols <- c(display_cols, "source")
        }

        # Only keep columns that exist
        display_cols <- display_cols[display_cols %in% colnames(all_enrichment)]

        if (length(display_cols) > 0) {
          all_enrichment <- all_enrichment[, display_cols, drop = FALSE]

          # Format numeric columns for better readability
          if ("p.value" %in% colnames(all_enrichment)) {
            all_enrichment$p.value <- formatC(all_enrichment$p.value, format = "e", digits = 3)
          }
          if ("p.adjust" %in% colnames(all_enrichment)) {
            all_enrichment$p.adjust <- formatC(all_enrichment$p.adjust, format = "e", digits = 3)
          }
          if ("RichFactor" %in% colnames(all_enrichment)) {
            all_enrichment$RichFactor <- round(all_enrichment$RichFactor, 4)
          }
          if ("FoldEnrichment" %in% colnames(all_enrichment)) {
            all_enrichment$FoldEnrichment <- round(all_enrichment$FoldEnrichment, 3)
          }
          if ("zscore" %in% colnames(all_enrichment)) {
            all_enrichment$zscore <- round(all_enrichment$zscore, 3)
          }
        }

        DT::datatable(all_enrichment,
                      options = list(pageLength = 25,
                                    scrollX = TRUE,
                                    scrollY = "500px",
                                    scrollCollapse = TRUE),
                      rownames = FALSE,
                      filter = "top",
                      caption = htmltools::tags$caption(
                        style = 'caption-side: top; text-align: left; color: #333; font-size: 14px;',
                        htmltools::strong('Enrichment Results: '),
                        'All columns from richR including KEGG pathway hierarchy (Level1/2/3) if applicable'
                      ))
      } else {
        DT::datatable(data.frame(Message = "No enrichment results found"))
      }
    } else {
      DT::datatable(data.frame(Message = "No enrichment results available"))
    }
  })

  # Download pathway results
  output$download_pathway <- downloadHandler(
    filename = function() {
      paste0("CellHet_Pathway_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
    },
    content = function(file) {
      req(values$pathway_results)

      withProgress(message = 'Preparing download...', value = 0, {
        incProgress(0.5, detail = "Exporting results")

        tryCatch({
          CellHet::exportEnrichmentResults(
            enrichment_results = values$pathway_results,
            file_path = file,
            format = "xlsx"
          )
          incProgress(1, detail = "Done!")
        }, error = function(e) {
          showNotification(
            paste("Error creating download:", e$message),
            type = "error",
            duration = 10
          )
        })
      })
    }
  )

  # ==================== PDF Download Handlers ====================

  # Download volcano plot as PDF
  output$download_volcano_pdf <- downloadHandler(
    filename = function() {
      paste0("CellHet_Volcano_", gsub("__", "_", input$selected_comparison), "_", format(Sys.time(), "%Y%m%d"), ".pdf")
    },
    content = function(file) {
      req(values$deg_results, input$selected_comparison)

      # Get settings
      width <- if (!is.null(input$volcano_width)) input$volcano_width else 8
      height <- if (!is.null(input$volcano_height)) input$volcano_height else 6

      # Recreate the plot
      deg_data <- values$deg_results$degs[[input$selected_comparison]]
      deg_data$neg_log10_p <- -log10(deg_data$p_val_adj)
      deg_data$color <- "Not Significant"
      deg_data$color[deg_data$significant & deg_data$avg_log2FC > 0] <- "Up-regulated"
      deg_data$color[deg_data$significant & deg_data$avg_log2FC < 0] <- "Down-regulated"

      deg_data$label <- ""
      if (sum(deg_data$significant) > 0) {
        sig_data <- deg_data[deg_data$significant, ]
        sig_data <- sig_data[order(abs(sig_data$avg_log2FC), decreasing = TRUE), ]
        top_genes <- head(sig_data$gene, 10)
        deg_data$label[deg_data$gene %in% top_genes] <- deg_data$gene[deg_data$gene %in% top_genes]
      }

      p_threshold <- if (!is.null(input$p_val_threshold) && input$pval_type == "unadjusted") {
        -log10(input$p_val_threshold)
      } else {
        -log10(input$p_val_adj_threshold)
      }

      font_size <- if (!is.null(input$volcano_fontsize)) input$volcano_fontsize else 12

      # Get color settings
      up_color <- if (!is.null(input$volcano_up_color)) input$volcano_up_color else "#E67E22"
      down_color <- if (!is.null(input$volcano_down_color)) input$volcano_down_color else "#00ACC1"
      ns_color <- if (!is.null(input$volcano_ns_color)) input$volcano_ns_color else "gray70"

      p <- ggplot(deg_data, aes(x = avg_log2FC, y = neg_log10_p, color = color)) +
        geom_point(alpha = 0.6, size = 2) +
        scale_color_manual(values = c("Up-regulated" = up_color,
                                     "Down-regulated" = down_color,
                                     "Not Significant" = ns_color)) +
        geom_hline(yintercept = p_threshold, linetype = "dashed", color = "red") +
        geom_vline(xintercept = c(-input$logfc_threshold, input$logfc_threshold),
                  linetype = "dashed", color = "red") +
        theme_minimal(base_size = font_size) +
        labs(title = paste("Volcano Plot:", input$selected_comparison),
             x = "Log2 Fold Change",
             y = "-Log10 Adjusted P-value",
             color = "Regulation")

      if (any(deg_data$label != "")) {
        p <- p + ggrepel::geom_text_repel(aes(label = label), size = font_size / 4, max.overlaps = 20)
      }

      pdf(file, width = width, height = height)
      print(p)
      dev.off()
    }
  )

  # Download heatmap as PDF
  output$download_heatmap_pdf <- downloadHandler(
    filename = function() {
      paste0("CellHet_Heatmap_", format(Sys.time(), "%Y%m%d"), ".pdf")
    },
    content = function(file) {
      req(values$deg_results)

      width <- if (!is.null(input$heatmap_width)) input$heatmap_width else 10
      height <- if (!is.null(input$heatmap_height)) input$heatmap_height else 8
      font_size <- if (!is.null(input$heatmap_fontsize)) input$heatmap_fontsize else 12

      # Get color settings
      low_col <- if (!is.null(input$heatmap_low_color)) input$heatmap_low_color else "#4DA8A0"
      mid_col <- if (!is.null(input$heatmap_mid_color)) input$heatmap_mid_color else "white"
      high_col <- if (!is.null(input$heatmap_high_color)) input$heatmap_high_color else "#E67E22"

      pdf(file, width = width, height = height)
      print(CellHet::visualizeDEGHeatmap(
        deg_results = values$deg_results,
        comparison_col = "comparison",
        cell_type_col = "cell_type",
        show_counts = TRUE,
        cluster_rows = TRUE,
        low_color = low_col,
        mid_color = mid_col,
        high_color = high_col,
        title = paste("Top", input$n_top_genes, "DEGs Across Cell Types and Comparisons"),
        title_size = font_size + 2,
        count_size = font_size / 4,
        top_label_size = font_size / 3
      ))
      dev.off()
    }
  )

  # Download upset plot as PDF
  output$download_upset_pdf <- downloadHandler(
    filename = function() {
      paste0("CellHet_UpSet_", format(Sys.time(), "%Y%m%d"), ".pdf")
    },
    content = function(file) {
      req(values$deg_results)

      width <- if (!is.null(input$upset_width)) input$upset_width else 10
      height <- if (!is.null(input$upset_height)) input$upset_height else 6
      font_size <- if (!is.null(input$upset_fontsize)) input$upset_fontsize else 12

      # Get color settings
      bar_color <- if (!is.null(input$upset_bar_color)) input$upset_bar_color else "#2670B8"
      point_color <- if (!is.null(input$upset_point_color)) input$upset_point_color else "black"
      use_palette <- if (!is.null(input$upset_use_palette)) input$upset_use_palette else TRUE

      # Set sets_bar_colors based on checkbox
      sets_colors <- if (use_palette) NULL else bar_color

      pdf(file, width = width, height = height)
      print(CellHet::visualizeDEGUpset(
        deg_results = values$deg_results,
        by_cell_type = TRUE,
        direction = "both",
        min_size = 3,
        bar_width = 0.6,
        main_bar_color = bar_color,
        intersection_color = point_color,
        point_outline_color = point_color,
        sets_bar_colors = sets_colors,
        title = "DEG Overlap Across Cell Types",
        text_size = font_size,
        set_text_size = font_size,
        set_label_size = font_size / 4,
        intersection_label_size = font_size / 4
      ))
      dev.off()
    }
  )

  # Download pathway barplot as PDF
  output$download_pathway_bar_pdf <- downloadHandler(
    filename = function() {
      paste0("CellHet_Pathway_Bar_", format(Sys.time(), "%Y%m%d"), ".pdf")
    },
    content = function(file) {
      req(values$pathway_results)

      width <- if (!is.null(input$pathway_bar_width)) input$pathway_bar_width else 10
      height <- if (!is.null(input$pathway_bar_height)) input$pathway_bar_height else 8

      if (!is.null(values$pathway_results$plots) && "bar" %in% names(values$pathway_results$plots)) {
        pdf(file, width = width, height = height)
        print(values$pathway_results$plots$bar)
        dev.off()
      }
    }
  )

  # Download pathway dotplot as PDF
  output$download_pathway_dot_pdf <- downloadHandler(
    filename = function() {
      paste0("CellHet_Pathway_Dot_", format(Sys.time(), "%Y%m%d"), ".pdf")
    },
    content = function(file) {
      req(values$pathway_results)

      width <- if (!is.null(input$pathway_dot_width)) input$pathway_dot_width else 10
      height <- if (!is.null(input$pathway_dot_height)) input$pathway_dot_height else 8

      if (!is.null(values$pathway_results$plots) && "dot" %in% names(values$pathway_results$plots)) {
        pdf(file, width = width, height = height)
        print(values$pathway_results$plots$dot)
        dev.off()
      }
    }
  )

  # Download pathway network as PDF
  output$download_pathway_net_pdf <- downloadHandler(
    filename = function() {
      paste0("CellHet_Pathway_Network_", format(Sys.time(), "%Y%m%d"), ".pdf")
    },
    content = function(file) {
      req(values$pathway_results)

      width <- if (!is.null(input$pathway_net_width)) input$pathway_net_width else 12
      height <- if (!is.null(input$pathway_net_height)) input$pathway_net_height else 10
      font_size <- if (!is.null(input$pathway_net_fontsize)) input$pathway_net_fontsize else 12
      node_color <- if (!is.null(input$pathway_net_node_color)) input$pathway_net_node_color else "#2670B8"

      if (!is.null(values$pathway_results$plots) && "network" %in% names(values$pathway_results$plots)) {
        pdf(file, width = width, height = height)
        # Display richR network plot as-is (color modifications not supported)
        print(values$pathway_results$plots$network)
        dev.off()
      }
    }
  )

  # ==================== Trajectory Analysis ====================

  # Trajectory variable selector
  output$trajectory_var_ui <- renderUI({
    req(values$metadata)

    # Find numeric columns for trajectory
    numeric_cols <- names(values$metadata)[sapply(values$metadata, is.numeric)]

    if (length(numeric_cols) == 0) {
      helpText("No numeric columns found for trajectory analysis")
    } else {
      selectInput("trajectory_var", "Trajectory Variable:",
                 choices = numeric_cols,
                 selected = if("pseudotime" %in% numeric_cols) "pseudotime" else numeric_cols[1])
    }
  })

  # Run trajectory analysis
  observeEvent(input$run_trajectory, {
    req(values$seurat_obj, input$trajectory_var)

    withProgress(message = 'Running trajectory analysis...', value = 0, {

      incProgress(0.3, detail = "Analyzing trajectory")

      tryCatch({
        values$trajectory_results <- CellHet::analyzeDEGTrajectory(
          object = values$seurat_obj,
          trajectory_var = input$trajectory_var,
          cell_type_var = input$cell_type_var,
          group_var = NULL,
          n_bins = input$n_bins,
          min_cells_per_bin = input$min_cells_per_bin,
          smooth_method = input$smooth_method,
          n_top_genes = input$n_top_trajectory_genes,
          cluster_genes = TRUE,
          plot_type = "both"
        )

        incProgress(1, detail = "Done!")
        showNotification("Trajectory analysis completed!", type = "message")

      }, error = function(e) {
        showNotification(paste("Error in trajectory analysis:", e$message), type = "error", duration = 10)
      })
    })
  })

  # Trajectory heatmap
  output$trajectory_heatmap <- renderPlot({
    req(values$trajectory_results)

    if (!is.null(values$trajectory_results$plots$heatmap)) {
      # Get color settings
      low_col <- if (!is.null(input$traj_low_color)) input$traj_low_color else "deepskyblue"
      mid_col <- if (!is.null(input$traj_mid_color)) input$traj_mid_color else "white"
      high_col <- if (!is.null(input$traj_high_color)) input$traj_high_color else "darkorange"

      # Check if it's a ggplot object (fallback heatmap) or ComplexHeatmap
      if (inherits(values$trajectory_results$plots$heatmap, "gg")) {
        # It's a ggplot, we can modify the color scale
        values$trajectory_results$plots$heatmap +
          scale_fill_gradient2(low = low_col, mid = mid_col, high = high_col,
                              na.value = "grey90")
      } else {
        # It's a ComplexHeatmap, display as-is (colors are baked in)
        # Note: ComplexHeatmap objects need special handling to display in Shiny
        grid::grid.newpage()
        ComplexHeatmap::draw(values$trajectory_results$plots$heatmap)
      }
    } else {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                label = "Heatmap not available",
                size = 6) +
        theme_void()
    }
  })

  # Gene selector for line plots
  output$gene_select_ui <- renderUI({
    req(values$trajectory_results)

    top_genes <- values$trajectory_results$top_genes
    selectInput("selected_genes", "Select Genes to Plot:",
               choices = top_genes,
               selected = head(top_genes, 5),
               multiple = TRUE)
  })

  # Trajectory line plot
  output$trajectory_lineplot <- renderPlot({
    req(values$trajectory_results, input$selected_genes)

    tryCatch({
      CellHet::plotGeneTrajectory(
        trajectory_results = values$trajectory_results,
        genes = input$selected_genes,
        plot_type = "line",
        add_loess_fit = FALSE
      )
    }, error = function(e) {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                label = paste("Error:", e$message),
                size = 5) +
        theme_void()
    })
  })

  # Trajectory summary
  output$trajectory_summary <- renderPrint({
    req(values$trajectory_results)

    cat("Trajectory Analysis Summary\n")
    cat("===========================\n\n")
    cat("Trajectory Variable:", values$trajectory_results$trajectory_info$trajectory_var, "\n")
    cat("Trajectory Range:",
        round(values$trajectory_results$trajectory_info$min_value, 3), "to",
        round(values$trajectory_results$trajectory_info$max_value, 3), "\n")
    cat("Number of Bins:", length(values$trajectory_results$trajectory_info$bin_counts), "\n")
    cat("Cells per Bin:", paste(values$trajectory_results$trajectory_info$bin_counts, collapse = ", "), "\n\n")
    cat("Top Variable Genes:", length(values$trajectory_results$top_genes), "\n")
    cat("Gene Clusters:",
        if(!is.null(values$trajectory_results$gene_clusters)) length(unique(values$trajectory_results$gene_clusters)) else "None",
        "\n")
  })

  # Trajectory bins plot
  output$trajectory_bins_plot <- renderPlot({
    req(values$trajectory_results)

    bin_data <- data.frame(
      Bin = 1:length(values$trajectory_results$trajectory_info$bin_counts),
      Count = as.vector(values$trajectory_results$trajectory_info$bin_counts)
    )

    ggplot(bin_data, aes(x = factor(Bin), y = Count)) +
      geom_bar(stat = "identity", fill = "#2670B8") +
      theme_minimal() +
      labs(x = "Bin", y = "Number of Cells", title = "Cells per Trajectory Bin")
  })

  # ==================== Heterogeneity Analysis ====================

  # Run heterogeneity analysis
  observeEvent(input$run_het, {
    req(values$deg_results)

    # Check if pathway results are needed
    if (input$het_level %in% c("pathway", "both") && is.null(values$pathway_results)) {
      showNotification("Pathway results required for pathway-level heterogeneity. Please run pathway analysis first.",
                      type = "warning", duration = 10)
      return(NULL)
    }

    withProgress(message = 'Calculating heterogeneity...', value = 0, {

      incProgress(0.3, detail = "Computing similarity matrix")

      tryCatch({
        values$het_results <- CellHet::visualizeHeterogeneity(
          deg_results = values$deg_results,
          pathway_results = values$pathway_results,
          layout_type = input$het_layout,
          level = input$het_level,
          similarity_metric = input$het_metric
        )

        incProgress(1, detail = "Done!")
        showNotification("Heterogeneity analysis completed!", type = "message")

      }, error = function(e) {
        showNotification(paste("Error in heterogeneity analysis:", e$message), type = "error", duration = 10)
      })
    })
  })

  # Heterogeneity plot
  output$het_plot <- renderPlot({
    req(values$het_results)

    if (!is.null(values$het_results$plot)) {
      # Check if it's a ComplexHeatmap or ggplot object
      if (inherits(values$het_results$plot, "Heatmap")) {
        # It's a ComplexHeatmap object
        grid::grid.newpage()
        ComplexHeatmap::draw(values$het_results$plot)
      } else if (inherits(values$het_results$plot, "gg")) {
        # It's a ggplot object
        values$het_results$plot
      } else {
        # For other plot types (like igraph), just print it
        print(values$het_results$plot)
      }
    } else {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                label = "No plot available. Please run heterogeneity analysis first.",
                size = 6) +
        theme_void()
    }
  })

  # Heterogeneity matrix table
  output$het_matrix_table <- DT::renderDataTable({
    req(values$het_results)

    if (!is.null(values$het_results$heterogeneity_matrix)) {
      # Convert matrix to data frame with row names as a column
      het_df <- as.data.frame(values$het_results$heterogeneity_matrix)
      het_df <- cbind(CellType = rownames(het_df), het_df)
      rownames(het_df) <- NULL

      DT::datatable(het_df,
                    options = list(pageLength = 10, scrollX = TRUE, digits = 3),
                    rownames = FALSE) %>%
        DT::formatRound(columns = 2:ncol(het_df), digits = 3)
    } else {
      DT::datatable(data.frame(Message = "No heterogeneity matrix available"))
    }
  })

  # Similarity matrix table
  output$het_similarity_table <- DT::renderDataTable({
    req(values$het_results)

    if (!is.null(values$het_results$similarity_matrix)) {
      # Convert matrix to data frame with row names as a column
      sim_df <- as.data.frame(values$het_results$similarity_matrix)
      sim_df <- cbind(CellType = rownames(sim_df), sim_df)
      rownames(sim_df) <- NULL

      DT::datatable(sim_df,
                    options = list(pageLength = 10, scrollX = TRUE, digits = 3),
                    rownames = FALSE) %>%
        DT::formatRound(columns = 2:ncol(sim_df), digits = 3)
    } else {
      DT::datatable(data.frame(Message = "No similarity matrix available"))
    }
  })

  # Session info for Help tab
  output$session_info <- renderPrint({
    sessionInfo()
  })
}

# Run the application
shinyApp(ui = ui, server = server)
