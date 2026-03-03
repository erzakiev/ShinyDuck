mod_pca_ui <- function(id) {
  ns <- NS(id)

  tagList(
    tabsetPanel(
      tabPanel("Top 500 HVG",
               plotOutput(ns("pca_top500"), height = 800)
      ),
      tabPanel("All genes",
               plotOutput(ns("pca_all"), height = 800)
      )
    )
  )
}

mod_pca_server <- function(id, rv) {

  moduleServer(id, function(input, output, session) {
    vst_data <- reactive({
      req(rv$app_state == 'ready')
      DESeq2::vst(rv$txi_deseq)
    })

    pca_all <- reactive({
      dsq <- vst_data()

      p <- DESeq2::plotPCA(dsq,
                   intgroup = "Group",
                   ntop = nrow(SummarizedExperiment::assay(dsq)))

      p + ggrepel::geom_label_repel(
        label = colnames(dsq),
        max.overlaps = 30
      ) +
        ggplot2::theme_classic()
    })

    pca_top500 <- reactive({
      dsq <- vst_data()

      p <- DESeq2::plotPCA(dsq,
                   intgroup = "Group",
                   ntop = 500)

      p + ggrepel::geom_label_repel(
        label = colnames(dsq),
        max.overlaps = 30
      ) +
        ggplot2::theme_classic()
    })

    output$pca_all <- renderPlot({
      pca_all()
    })

    output$pca_top500 <- renderPlot({
      pca_top500()
    })

  })
}
