#' Run the RhodopXin Shiny App
#'
#' @description This function runs the Shiny web application for this package,
#'   providing an interactive interface for a common workflow provided by
#'   RhodopXIn. This function does not return, but a webpage will be opened.
#'
#'
#' @examples
#' \dontrun{
#' runRhodopXin()
#' }
#'
#' @export
#' @importFrom shiny runApp
runRhodopXin <- function() {
  # Run the Shiny application from its directory
  app_dir <- system.file("shiny-scripts", package = "RhodopXin")
  shiny::runApp(appDir = app_dir, display.mode = "normal")

  # Nothing to return
  return(invisible(NULL))
}

# [END]
