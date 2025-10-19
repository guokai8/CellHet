#' Launch CellHet Shiny Application
#'
#' This function launches the interactive Shiny application for CellHet analysis.
#' The app provides a user-friendly interface for performing differential expression,
#' pathway enrichment, trajectory analysis, and heterogeneity quantification.
#'
#' @param launch.browser Logical, whether to launch the app in a browser (default: TRUE)
#' @param port Integer, the port to run the app on (default: NULL, uses random port)
#' @param host Character, the host to run the app on (default: "127.0.0.1")
#'
#' @return Launches the Shiny application
#'
#' @examples
#' \dontrun{
#' # Launch the app
#' runCellHetApp()
#'
#' # Launch on specific port
#' runCellHetApp(port = 3838)
#' }
#'
#' @export
runCellHetApp <- function(launch.browser = TRUE, port = NULL, host = "127.0.0.1") {

  # Check if required packages are installed
  required_packages <- c("shiny", "shinythemes", "DT", "ggplot2")

  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop(paste0(
      "The following packages are required to run the Shiny app:\n",
      paste(missing_packages, collapse = ", "),
      "\n\nPlease install them using:\n",
      "install.packages(c(", paste0("'", missing_packages, "'", collapse = ", "), "))"
    ))
  }

  # Get app directory
  app_dir <- system.file("shiny", package = "CellHet")

  if (app_dir == "") {
    stop("Could not find Shiny app directory. Please ensure the package is properly installed.")
  }

  # Launch options
  launch_options <- list(
    launch.browser = launch.browser,
    host = host
  )

  if (!is.null(port)) {
    launch_options$port <- port
  }

  # Launch the app
  message("Launching CellHet Shiny App...")
  message("Use Ctrl+C (or Cmd+C on Mac) to stop the application")

  shiny::runApp(
    appDir = app_dir,
    launch.browser = launch_options$launch.browser,
    host = launch_options$host,
    port = if (!is.null(port)) port else NULL
  )
}
