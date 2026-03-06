mod_deseq_deg_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(ns("text_DEGs"))),
    h4(textOutput(ns("text_which_is_Control_in_DEG"))),
    DT::DTOutput(ns("DESeq_DEGs")),
    uiOutput(ns("DESeq_DEGsMultitab"))
  )
}

mod_deseq_deg_server <- function(id, rv) {

  moduleServer(id, function(input, output, session) {

    output$DESeq_DEGs <- DT::renderDT({
      req(rv$app_state == 'ready')
      req(rv$multiple_groups == 0)
      message("Preparing DEG table for mono-group")
      rstxdsq <- rv$res_txi_deseq
      DT::datatable(as.data.frame(rstxdsq), filter = 'top', extensions = 'Buttons',
                options = list(
                  dom = "Bl<'search'>rtip",
                  buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                  searchCols = list(NULL, NULL, NULL, NULL,
                                    NULL, NULL, NULL, NULL, NULL, NULL, NULL, list(search = 'Significant (padj < 0.05)'))
                )#,
                #callback = JS(callback)
                )
    },
    server=FALSE, rownames = FALSE
    )

    output$DESeq_DEGsMultitab <- renderUI({
      req(rv$app_state == 'ready')
      req(rv$multiple_groups == 1)

      message("Preparing DEG table for multi-groups")

      nTabs = length(rv$res_DEGs_txi_deseq)
      rstxdsq <- rv$res_txi_deseq

      myTabs = lapply(1: nTabs, function(x){
        tabPanel(paste(strsplit(names(rstxdsq[x]), split = "_")[[1]][2:4], collapse=' '),
                 DT::renderDT({DT::datatable(rstxdsq[[x]], filter = 'top', extensions = 'Buttons',
                                     options = list(dom = "Blrtip",
                                                    buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                                    searchCols = list(NULL, NULL, NULL, NULL,
                                                                      NULL, NULL, NULL, NULL, NULL, NULL, NULL, list(search = 'Significant (padj < 0.05)'))
                                     ), #callback = JS(callback),
                 )},
                 server=FALSE, rownames = FALSE))});
      return(do.call(tabsetPanel, myTabs))
    })

    output$text_DEGs <- renderText({
      req(rv$app_state == 'ready')
      "List of differentially expressed genes:"
    })
    output$text_which_is_Control_in_DEG <- renderText({
      req(rv$app_state == 'ready')
      paste0("The control (background) is ", levels(rv$txi_deseq_deseq$Group)[1])
    })

    n_contrasts <- reactive({
      req(rv$app_state == 'ready')
      return(length(resultsNames(rv$res_DEGs_txi_deseq)))
    })

  })
}
