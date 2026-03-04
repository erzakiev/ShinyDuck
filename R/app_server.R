
app_server <- function(input, output, session, data_root) {

  rv <- make_state()

  roots <- c(
    home = data_root
  )

  house_path <- data_root

  mod_project_setup_server(
    id = "setup",
    rv = rv,
    roots = roots,
    house_path = house_path
  )
  mod_fastqc_server("fastqc", rv)
  mod_tpm_table_server("tpm", rv)
  mod_deseq_deg_server("deg", rv)
  mod_pca_server("pca", rv)
  mod_volcano_server("volcano", rv)
  mod_coldata_server("coldata", rv)
  mod_enrich_server("enrich", rv)
  #mod_gsva_server("gsva", rv)
  #mod_download_server("download", rv)

  observeEvent(input$sidebar, {
    message("Active tab changed to: ", input$sidebar)
  }, ignoreNULL = FALSE)

  observeEvent(rv$app_state, {
    if(rv$app_state == "idle") {
      session$sendCustomMessage("toggle_menu_items", list(disable = TRUE))
    } else {
      session$sendCustomMessage("toggle_menu_items", list(disable = FALSE))
    }
  })

  observeEvent(rv$projFolderFull, {
    header <- basename(rv$projFolderFull)
    shinyjs::html("pageHeader", header)
  })
}
