mod_pca_ui <- function(id) {
  ns <- NS(id)

  tagList(
    tabsetPanel(
      tabPanel("Top 500 HVG",
               plotOutput(ns("pca_top500"), height = 800)#,
               #textInput('nGenes', placeholder = 'how many top contributors to PCs to show?'),
               #textOutput(ns('pca_top500_pc1')),
               #textOutput(ns('pca_top500_pc2')),
      ),
      tabPanel("All genes",
               plotOutput(ns("pca_all"), height = 800)#,
               #textInput('nGenes', placeholder = 'how many top contributors to PCs to show?'),
               #textOutput(ns('pca_all_pc1')),
               #textOutput(ns('pca_all_pc2')),
      )
    )
  )
}

mod_pca_server <- function(id, rv) {

  moduleServer(id, function(input, output, session) {
    pca_all <- reactive({
      dsq <- rv$vst_data

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
      dsq <- rv$vst_data

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


    #pca_500 <- reactive({prcomp(t(SummarizedExperiment::assay(rv$vst_data)), center = TRUE, scale. = FALSE)})

    #output$pca_top500_pc1 <- renderText({
    #  sort(abs(pca_500()$rotation[, 1]), decreasing = TRUE)[1:input$nGenes]
    #})

    #output$pca_top500_pc2 <- renderText({
    #  sort(abs(pca_500()$rotation[, 2]), decreasing = TRUE)[1:input$nGenes]
    #})
  })
}
