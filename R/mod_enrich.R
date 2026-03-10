mod_enrich_ui <- function(id) {
  ns <- NS(id)
  textOutput(ns("Enrichments_title"))
  tabsetPanel(id=ns('Enrichments_panel'), type = "tabs",
              tabPanel(ns("Dotplots"),
                       textInput(ns('showCategory'), label = 'Show how many top enriched terms?', value = "50", width='15%'),
                       uiOutput(ns('GO_dotplot_multitab')), plotOutput(ns("GO_dotplot"), height = '2000px')),
              tabPanel(ns("Tables of terms"), uiOutput(ns('GO_multitab')),  DT::DTOutput(ns("GO")))
  )
}

mod_enrich_server <- function(id, rv) {

  moduleServer(id, function(input, output, session) {

    output$Enrichments_title<- renderText({
      req(rv$app_state == 'ready')
      return("Enriched terms:")
    })

    output$GO<- DT::renderDT({
      req(rv$app_state == 'ready')
      req(rv$multiple_groups == 0)
      varbl <- format_df_numbers(as.data.frame(rv$GO_result))
      varbl$geneID <- gsub(pattern='/', replacement=' ', varbl$geneID)
      varbl
    }, rownames = FALSE,
    extensions = 'Buttons',
    options = list(
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
      scrollY = TRUE,
      scrollX = TRUE,
      deferRender = TRUE,
      scroller = TRUE,
      fixedHeader = TRUE),
    server=FALSE)


    output$GO_multitab <- renderUI({
      req(rv$app_state == 'ready')
      req(rv$multiple_groups == 1)

      nTabs = length(rv$GO_result)

      toOutput <- list()
      for (o in 1:nTabs){
        toOutput[[o]] <- format_df_numbers(as.data.frame(rv$GO_result[[o]]))
        rownames(toOutput[[o]]) <- NULL
        toOutput[[o]]$geneID <- gsub(pattern='/', replacement=' ', toOutput[[o]]$geneID)
      }

      myTabs = lapply(1: nTabs, function(x){tabPanel(names(rv$res_txi_deseq)[x] |> gsub(pattern = 'Group_', replacement = ''),
                                                    DT::renderDT(DT::datatable(toOutput[[x]], filter = 'top', rownames = FALSE,extensions = 'Buttons',
                                                                          options = list(dom = "Blrtip",
                                                                                         buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                                                                         scrollY = TRUE,
                                                                                         scrollX = TRUE,
                                                                                         deferRender = TRUE,
                                                                                         scroller = TRUE,
                                                                                         fixedHeader = TRUE)),
                                                                server=FALSE))});
      return(do.call(tabsetPanel, myTabs))
    })

    output$GO_dotplot<- renderPlot({
      req(rv$app_state == 'ready')
      req(rv$multiple_groups == 0)
      ncat <- as.numeric(input$showCategory)

      p <- enrichplot::dotplot(rv$GO_result,
                               showCategory=50,
                               x='p.adjust',
                               decreasing=F)
      png(file.path(rv$projFolderFull,'GO_dotplot_300dpi.png'),
          width = 10,
          height = 20,
          res = 300, units = 'in')
      print(p)
      dev.off()

      if(ncat!=50){
        p_custom <- enrichplot::dotplot(rv$GO_result,
                                        showCategory=ncat,
                                        x='p.adjust',
                                        decreasing=F)
        png(file.path(rv$projFolderFull,'GO_dotplot_custom_number_of_Terms_300dpi.png'),
            width = 10,
            height = round(20*ncat/50, 1),
            res = 300, units = 'in')
        print(p_custom)
        dev.off()
        click("downloadDotPlot")
        return(p_custom)
      }
      return(p)
    })

    output$GO_dotplot_multitab <- renderUI({
      req(rv$app_state == 'ready')
      req(rv$multiple_groups == 1)
      p <- list()
      tytl <- list()
      ncat <- as.numeric(input$showCategory)
      for (i in 1:length(rv$res_txi_deseq)){
        tytl[[i]] <- names(rv$res_txi_deseq)[i] |> gsub(pattern = 'Group_', replacement = '')
        if(is.null(rv$GO_result[[i]])){ p[[i]] <- ggplot()+theme_void()} else {
          if(ncat!=50){
            p[[i]] <- enrichplot::dotplot(rv$GO_result[[i]],
                                          showCategory=ncat,
                                          x='p.adjust',
                                          decreasing=F)
            flnm <- file.path(rv$projFolderFull,paste0('GO_dotplot_custom_number_of_Terms_300dpi_', tytl[[i]],'.png'))
            png(filename = flnm,
                width = 10,
                height = round(20*ncat/50, 1),
                res = 300, units = 'in')
            print(p[[i]])
            dev.off()
          } else {
            p[[i]] <- enrichplot::dotplot(rv$GO_result[[i]],
                                          showCategory=50,
                                          x='p.adjust',
                                          decreasing=F)
            flnm <- file.path(rv$projFolderFull,paste0('GO_dotplot_300dpi_', tytl[[i]],'.png'))
            png(filename = flnm,
                width = 10,
                height = 20,
                res = 300, units = 'in')
            print(p[[i]])
            dev.off()
          }

        }
      }
      nTabs = length(rv$res_txi_deseq)
      myTabs = lapply(1: nTabs, function(x){tabPanel(tytl[[x]],
                                                     renderPlot(p[[x]], height = 2000))});

      return(do.call(tabsetPanel, myTabs))
    })

  })
}
