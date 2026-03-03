mod_volcano_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(ns("Volcano_title"))),
    plotly::plotlyOutput(ns("Volcano"), height = 900)
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

      df$SignificanceLevel[df$padj < 0.05 & abs(df$log2FoldChange) > 0.5] <- "Significant&FoldChange"
      df$SignificanceLevel[df$padj > 0.05 & abs(df$log2FoldChange) > 0.5] <- "FoldChange"
      df$SignificanceLevel[df$padj < 0.05 & abs(df$log2FoldChange) <= 0.5] <- "Significant"

      df$padj[df$padj == 0] <- min(df$padj[df$padj != 0]) / 10000

      df
    })

    volcano_plot <- reactive({
      df <- volcano_df()

      tx <- highlight_key(df, ~SYMBOL)

      base <- plot_ly(tx, height = 900) %>%
        add_trace(
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
        ) %>%
        add_trace(
          data = df %>%
            dplyr::filter(SignificanceLevel == "Significant&FoldChange") %>%
            dplyr::slice_min(padj, n = 20),
          x = ~log2FoldChange,
          y = ~-log10(padj),
          text = ~SYMBOL,
          mode = "text",
          textposition = "top right",
          name = "Annotations"
        ) %>%
        layout(
          xaxis = list(showgrid = FALSE),
          yaxis = list(showgrid = FALSE),
          title = ""
        )

      highlight(
        base,
        on = "plotly_click",
        selectize = TRUE,
        dynamic = TRUE,
        persistent = TRUE,
        opacityDim = 0.07
      )
    })

    output$Volcano <- plotly::renderPlotly({
      volcano_plot()
    })

  })
}
