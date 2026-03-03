#' Run BulkRNA Suite App
#' @export

run_app <- function(data_root = "/media/minicluster/Data/FASTQ") {

  shiny::shinyApp(
    ui = app_ui(),
    server = function(input, output, session) {
      app_server(input, output, session, data_root)
    }
  )
}
