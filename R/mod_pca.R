mod_pca_ui <- function(id) {
  ns <- NS(id)

  tagList(
    tabsetPanel(
      tabPanel("Top 500 HVG",
               column(9,
                      plotOutput(ns("pca_top500"), height = 800),
                      ),
               column(3,
                      h4('Top contributing genes'),
                      textInput(ns('nGenes'),
                                placeholder = '10',
                                value = 10,
                                label = 'how many top contributors to PCs to show?'),
                      DT::DTOutput(ns('pca_top500_genes'))
               )
      ),
      tabPanel("All genes",
               column(9,
                      plotOutput(ns("pca_all"), height = 800),
                      ),
               column(3,
                      h4('Top contributing genes'),
                      textInput(ns('nGenes'),
                                placeholder = '10',
                                value = 10,
                                label = 'how many top contributors to PCs to show?'),
                      DT::DTOutput(ns('pca_all_genes'))
               )
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
      p <- pca_all()
      png(file.path(rv$projFolderFull, 'pca_all.png'), width = 10, height = 10, units = 'in', res = 300)
      print(p)
      dev.off()
      return(p)
    })

    output$pca_top500 <- renderPlot({
      p <- pca_top500()
      png(file.path(rv$projFolderFull, 'pca_top500.png'), width = 10, height = 10, units = 'in', res = 300)
      print(p)
      dev.off()
      return(p)
    })


    vsd <- reactive({SummarizedExperiment::assay(rv$vst_data)})

    pca_all_prcomp <- reactive({prcomp(t(convertDfRownamesEns2Symb(vsd(),rv$OrgDeeBee)), center = TRUE, scale. = FALSE)})

    output$pca_all_genes <- DT::renderDT({

      PC1_genes <- names(sort(abs(pca_all_prcomp()$rotation[, 1]), decreasing = TRUE)[1:input$nGenes])
      PC2_genes <- names(sort(abs(pca_all_prcomp()$rotation[, 2]), decreasing = TRUE)[1:input$nGenes])

      DT::datatable(
        data = data.frame(PC1=PC1_genes, PC2=PC2_genes),
        extensions = "Buttons",
        options = list(
          dom = "Blfrtip",
          buttons = c("copy", "csv", "excel", "pdf", "print"),
          pageLength = 25,
          scrollX = F
        ),
        rownames = FALSE
      )

    }, server = FALSE)

    output$pca_top500_genes <- DT::renderDT({

      rwvrs <- MatrixGenerics::rowVars(vsd())

      # Select top 500 genes by variance (exactly what plotPCA does)
      select <- order(rwvrs, decreasing = TRUE)[1:500]

      # Subset the matrix to only these genes
      dood <- prcomp(t(convertDfRownamesEns2Symb(vsd(),rv$OrgDeeBee)[select, ]), center = TRUE, scale. = FALSE)

      PC1_genes <- names(sort(abs(dood$rotation[, 1]), decreasing = TRUE)[1:input$nGenes])
      PC2_genes <- names(sort(abs(dood$rotation[, 2]), decreasing = TRUE)[1:input$nGenes])

      DT::datatable(
        data = data.frame(PC1=PC1_genes, PC2=PC2_genes),
        extensions = "Buttons",
        options = list(
          dom = "Blfrtip",
          buttons = c("copy", "csv", "excel", "pdf", "print"),
          pageLength = 25,
          scrollX = F
        ),
        rownames = FALSE
      )

    }, server = FALSE)
  })
}
