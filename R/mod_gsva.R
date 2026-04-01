mod_gsva_ui <- function(id) {
  ns <- NS(id)

  tagList(
    h3("GSVA Analysis for current groups of samples"),
    fluidRow(
      column(8, textInput(ns('genes_for_gsva'),
                          label = '',
                          value = "",
                          placeholder = 'Enter/paste a list of genes in SYMBOL and/or ENSEMBL formats',
                          width='100%')),
      column(1, div( style = "margin-top: 20px;",
                     actionButton(inputId = ns("runGSVA"), "Calculate")))
    ),
    verbatimTextOutput(ns("GSVA_genes_matching")),
    plotOutput(ns("GSVAplot"), height = '1200px'),
    DT::DTOutput(ns("GSVAtable"))
  )
}

mod_gsva_server <- function(id, rv) {

  moduleServer(id, function(input, output, session) {

    gsva_string <- reactive({
      req(rv$app_state=='ready')
      req(input$runGSVA)
      if(nchar(input$genes_for_gsva) == 0) return(NULL)
      toRet <- strsplits(input$genes_for_gsva, c(" ", ",", "/"))
      return(toRet)
    })

    gsva_string_matched <- reactiveVal('')


    output$GSVA_genes_matching<- renderText({
      req(rv$app_state=='ready')
      req(input$runGSVA)
      if(is.null(gsva_string())) return("List of genes shouldn't be empty")


      if(length(gsva_string_matched())<2){
        return(paste0('Only 1 gene was correctly mapped: ', gsva_string_matched()))
      } else {
        nnmpdgns <- setdiff(gsva_string(),
                            gsva_string_matched())
        if(length(nnmpdgns)==0){
          return(paste0('All genes were mapped correctly! The genes were: ',
                        paste(gsva_string_matched(), collapse = ",")))
        } else {
          return(paste0('Correctly mapped genes: ', paste(gsva_string_matched(), collapse = ","), '\n',
                        'Non-mapped genes: ', paste(nnmpdgns, collapse = ","), '\n',
                        "(If the list of unmapped genes is too long, try to replace them with their ENSEMBL notation)"
          ))
        }
      }
    })

    gsva_result <- reactive({
      req(rv$app_state=='ready')
      req(input$runGSVA)
      if(is.null(gsva_string())) return({})

      abund <- rv$txi_tpms
      groups <- abund[1,3:ncol(abund)]
      abund <- abund[-1,]

      abund <- abund[which(!duplicated(abund$ENSEMBL)),]
      rownames(abund) <- abund$ENSEMBL
      matching_genes <- unique(gsva_string()[which(toupper(gsva_string()) %in% toupper(c(abund$SYMBOL, abund$ENSEMBL)))])
      non_matching_genes <- unique(gsva_string()[which(!toupper(gsva_string()) %in% toupper(c(abund$SYMBOL, abund$ENSEMBL)))])

      if(length(matching_genes) < 2){
        return(NULL)
      } else {

        gsva_string_matched(matching_genes)

        which(toupper(c(abund$SYMBOL, abund$ENSEMBL)) %in% toupper(gsva_string())) -> v1

        ifelse(test=v1>nrow(abund), yes=v1-nrow(abund), no=v1) -> v2

        signature <- rownames(abund)[unique(v2)]
        abund <- abund[,-2]
        abund <- abund[,-1]

        gsva <- GSVA::gsva(GSVA::gsvaParam(as.matrix(dplyr::mutate_all(abund, function(x) as.numeric(as.character(x)))),
                           list(userSignature=signature)))
        gsva <- rbind(groups, gsva)
        gsva <- t(gsva)
        colnames(gsva) <- c('Group','Value')
        gsva <- as.data.frame(gsva)
        gsva[,2] <- as.numeric(gsva[,2])

        openxlsx::write.xlsx(x = gsva, file = file.path(rv$projFolderFull,'gsva_table.xlsx'))
        #saveRDS(gsva, file = file.path(rv$projFolderFull, 'gsva_result.RDS'))
        return(gsva)
      }
    })


    output$GSVAplot<- renderPlot({
      req(rv$app_state=='ready')
      req(input$runGSVA)
      if(is.null(gsva_result())) return({})

      fntsize=12
      p1 <- ggplot2::ggplot(gsva_result(), ggplot2::aes(x=Group, y=Value, fill=Group)) +
        ggplot2::geom_boxplot(width=0.4)+
        ggplot2::ggtitle(paste0("GSVA analysis"), subtitle = "on TPM counts") +
        ggplot2::xlab("") +
        ggplot2::ylab("GSVA value")+
        ggplot2::theme_classic() +
        ggplot2::theme(legend.position = "none") + ggplot2::theme(axis.text = ggplot2::element_text(size = fntsize))+ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))

      flnm <- paste0("GSVAsignature.png")
      png(filename = file.path(rv$projFolderFull, 'GSVAsignature.png'), width = 8, height = 6, units = "in", res =300)
      print(p1)
      dev.off()




      return(p1)
    })

    output$GSVAtable<- DT::renderDT({
      req(rv$app_state=='ready')
      req(input$runGSVA)
      if(is.null(gsva_result())) return(NULL)
      as.data.frame(gsva_result())
    }, extensions = 'Buttons',
    options = list(
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
    server=FALSE, rownames=T)

  })
}
