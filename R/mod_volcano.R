mod_volcano_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3('Interactive volcano plots'),
    h4(textOutput(ns("text_which_is_Control_in_DEG"))),
    h5('You can select genes on the plot or search for one in the search bar'),
    uiOutput(ns("Volcano_ui")),
    uiOutput(ns('VolcanoMultitab'), height=900)
  )
}

mod_volcano_server <- function(id, rv) {

  moduleServer(id, function(input, output, session) {

    volcano_df <- reactive({
      req(rv$app_state == 'ready')
      req(rv$multiple_groups == 0)
      message("Preparing volcano plot data for mono-group")

      df <- rv$res_txi_deseq

      df$SYMBOL[is.na(df$SYMBOL)] <- df$ENSEMBL

      df$SignificanceLevel <- "NS"

      df$SignificanceLevel[df$padj < 0.05 & abs(df$log2FoldChange) > 0.5] <- "Signif&FC"
      df$SignificanceLevel[df$padj > 0.05 & abs(df$log2FoldChange) > 0.5] <- "FC"
      df$SignificanceLevel[df$padj < 0.05 & abs(df$log2FoldChange) <= 0.5] <- "Signif"

      df$padj[df$padj == 0] <- min(df$padj[df$padj != 0]) / 10000

      df
    })

    output$Volcano <- reactive({
      df <- volcano_df()

      tx <- plotly::highlight_key(df, ~SYMBOL)

      base <- plotly::plot_ly(tx, height = 900) |>
        plotly::add_trace(
          x = ~log2FoldChange,
          y = ~-log10(padj),
          text = ~SYMBOL,
          mode = "markers",
          color = ~SignificanceLevel,
          colors = c("#32a852","black","blue","red"),
          hovertemplate = paste(
            "<b>%{text}</b><br>",
            "-log10(FDR): %{y:.2f}<br>",
            "log2FC: %{x:.2f}"
          )
        ) |>
        plotly::add_trace(
          data = df |>
            dplyr::filter(SignificanceLevel == "Signif&FC") |>
            dplyr::slice_min(padj, n = 20),
          x = ~log2FoldChange,
          y = ~-log10(padj),
          text = ~SYMBOL,
          mode = "text",
          textposition = "top right",
          name = "Annotations"
        ) |>
        plotly::layout(
          xaxis = list(showgrid = FALSE),
          yaxis = list(showgrid = FALSE),
          title = ""
        )

      plotly::highlight(
        base,
        on = "plotly_click",
        selectize = TRUE,
        dynamic = TRUE,
        persistent = TRUE,
        opacityDim = 0.07
      )
    })

    output$Volcano_ui <- renderUI({
      if (!isTRUE(rv$app_state == "ready") || rv$multiple_groups != 0)
        return(NULL)
      plotly::plotlyOutput(ns("Volcano"), height = 900)
    })


    output$VolcanoMultitab<- renderUI({

      req(rv$app_state == 'ready')
      req(rv$multiple_groups == 1)
      p <- list()
      tytl <- list()
      message("Preparing volcano plots data for multi-groups")
      for (i in 1:length(rv$res_txi_deseq)){
        rv$res_txi_deseq[[i]] -> df
        df[which(is.na(df$SYMBOL)),'SYMBOL'] <- df[which(is.na(df$SYMBOL)),'ENSEMBL']
        df$SignificanceLevel <- 'NS'
        df[which(df$padj < 0.05 & abs(df$log2FoldChange) > 0.5 ),"SignificanceLevel"] <- "Signif&FC"
        df[which(df$padj > 0.05 & abs(df$log2FoldChange) > 0.5 ),"SignificanceLevel"] <- "FC"
        df[which(df$padj < 0.05 & abs(df$log2FoldChange) < 0.5 ),"SignificanceLevel"] <- "Signif"

        tytl[[i]] <- names(rv$res_txi_deseq)[i] |> gsub(pattern = 'Group_', replacement = '')


        df.df <- as.data.frame(df)
        df.df$padj[which(df.df$padj==0)] <- min(df.df$padj[which(df.df$padj!=0)])/10000
        tx <- plotly::highlight_key(df.df, ~SYMBOL)

        # initiate a plotly object
        base <- plotly::plot_ly(tx, height=1000) |>
          plotly::add_trace(x = ~log2FoldChange,
                    y = ~-log10(padj),
                    text = ~SYMBOL, mode = 'markers',
                    color = ~SignificanceLevel,
                    colors = c("#32a852","black", "blue", "red"),
                    hovertemplate = paste('<b>%{text}</b><br>', '-log10(FDR): %{y:.2f}<br>','log2FC: %{x:.2f}'),
                    showlegend = T) |>
          plotly::add_trace(data = df.df |>
                      dplyr::filter(SignificanceLevel=='Signif&FC') |>
                        dplyr::top_n(-20, wt=padj),
                    x = ~log2FoldChange,
                    y = ~-log10(padj),
                    text = ~SYMBOL, mode = 'text',  textposition = "topright",
                    showlegend = T, name = 'Annotations') |>
          plotly::layout(xaxis=list(showgrid=F), yaxis=list(showgrid=F), title=tytl)

        p[[i]] <- plotly::highlight(
          base,
          on = "plotly_click",
          selectize = TRUE,
          dynamic = TRUE,
          persistent = TRUE,
          opacityDim = 0.07
        )
      }
      message(paste0('finished calculation of a volcano plot p[[i]] where i is ', i))



      nTabs = length(rv$res_txi_deseq)
      myTabs = lapply(1: nTabs, function(x){tabPanel(tytl[[x]], plotly::renderPlotly(p[[x]]))});

      return(do.call(tabsetPanel, myTabs))
    })

    output$text_which_is_Control_in_DEG <- renderText({
      req(rv$app_state == 'ready')
      paste0("The control (background) is ", levels(rv$txi_deseq_deseq$Group)[1])
    })

  })
}
