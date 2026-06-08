mod_dorothea_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Dorothea TF enrichment"),
    h4(''),
    plotOutput(ns("DorotheaPlot"), height = '1200px'),
    DT::DTOutput(ns("DorotheaTable"))
  )
}

mod_dorothea_server <- function(id, rv) {

  moduleServer(id, function(input, output, session) {
    #expr_mat <- assay(rv$vst_data)
    #library(dorothea)
    #library(dplyr)

#    data(dorothea_hs, package = "dorothea")

  })
}
