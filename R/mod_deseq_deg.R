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
      rstxdsq <- rv$res_txi_deseq |> dplyr::select(-c('lfcSE', 'stat', 'pvalue'))
      DT::datatable(format_df_numbers(as.data.frame(rstxdsq)), filter = 'top', extensions = 'Buttons',
                options = list(
                  dom = "Bl<'search'>rtip",
                  buttons = list(list(extend='copy'),
                                 list(extend='csv', filename = "DEG_table"),
                                 list(extend='excel', filename = "DEG_table"),
                                 list(extend='pdf', filename = "DEG_table"),
                                 list(extend='print'))
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
      rstxdsq <- lapply(rstxdsq, function(x) {
        dplyr::select(x, -c(lfcSE, stat, pvalue))
      })

      names(rstxdsq) <- names(rstxdsq) |> gsub(pattern = 'Group_', replacement = '')

      rstxdsq_altogether <- rstxdsq[[1]] |> dplyr::select(c("ENSEMBL", 'SYMBOL', 'baseMean'))
      for (counter in 1:length(rstxdsq)){
        rstxdsq_altogether <- cbind(rstxdsq_altogether,
                                    rstxdsq[[counter]] |> dplyr::select(c("log2FoldChange",
                                                                          'abs_log2FoldChange',
                                                                          'padj',
                                                                          'log10padj',
                                                                          'Significance')) |>
                                      dplyr::rename_with(~ paste0(.x, '_', names(rstxdsq)[counter])))
      }
      rstxdsq_altogether <- cbind(rstxdsq_altogether[,1:3],
                                  rstxdsq_altogether[, grep("^log2FoldChange", colnames(rstxdsq_altogether))],
                                  rstxdsq_altogether[, grep("^abs_log2FoldChange", colnames(rstxdsq_altogether))],
                                  rstxdsq_altogether[, grep("^padj", colnames(rstxdsq_altogether))],
                                  rstxdsq_altogether[, grep("^log10padj", colnames(rstxdsq_altogether))],
                                  rstxdsq_altogether[, grep("^Significance", colnames(rstxdsq_altogether))])


      rstxdsq[[length(rstxdsq)+1]] <- rstxdsq_altogether
      names(rstxdsq)[length(rstxdsq)] <- 'Altogether'

      #saveRDS(format_df_numbers(rstxdsq[[1]]), file = file.path(rv$projFolderFull, 'format_df_numbers(rstxdsq[[1]]).RDS'))

      myTabs = lapply(1: (nTabs+1), function(x){
        tabPanel(names(rstxdsq[x]),
                 DT::renderDT({DT::datatable(as.data.frame(format_df_numbers(rstxdsq[[x]])), filter = 'top', extensions = 'Buttons',
                                     options = list(dom = "Blrtip",
                                                    buttons = list(list(extend='copy'),
                                                                list(extend='csv', filename = paste0("DEG_table_", names(rstxdsq[x]))),
                                                                list(extend='excel', filename = paste0("DEG_table_", names(rstxdsq[x]))),
                                                                list(extend='pdf', filename = paste0("DEG_table_", names(rstxdsq[x]))),
                                                                list(extend='print')))
                                      #callback = JS(callback),
                 ) |>
                     color_numeric_column(data=as.data.frame(format_df_numbers(rstxdsq[[x]])),
                                        column="log2FoldChange")
                   },
                 server=FALSE, rownames = FALSE)

                 )});
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
