mod_project_setup_ui <- function(id) {
  ns <- NS(id)
  tabsetPanel(
    tabPanel(
      "Load previous analysis",
      fluidRow(
        column(3,
          shinyFiles::shinyDirButton(ns("infolder"), "Load previous project?", title="Select a folder")
        ),
        column(9, verbatimTextOutput(ns("previousprojname")))
      )
    ),
    tabPanel(
      "Run new analysis",
      fluidRow(
        column(3,
               textInput(ns("ProjectName"), "A new project name?", value = "..."),
               shinyFiles::shinyFilesButton(ns("infiles"), "Select .fastq files", title="Select fastq/fastq.gz", multiple=TRUE)
        ),
        column(3,
               shinyWidgets::pickerInput(
                 ns("referenceGenomeChoice"),
                 label = "Reference genome",
                 choices = c("Human GRCh38"=1, "Mouse GRCm39"=2),
                 selected = 2
               )
        )
      )
    )
  )
}

mod_project_setup_server <- function(id, rv, roots, house_path) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    shinyFiles::shinyFileChoose(input, "infiles", roots = roots, filetypes = c("", "gz","fastq","fq"))
    shinyFiles::shinyDirChoose(input, "infolder", roots = roots)


    observeEvent(input$infiles, {
      req(input$infiles)
      rv$files <- shinyFiles::parseFilePaths(roots, input$infiles)$datapath
      # rv$projFolderFull <- ... compute from files + ProjectName (your logic)
    })


    observeEvent(input$infolder, {

      req(input$infolder)

      inFolder <- paste0(
        house_path,
        paste0(input$infolder[[1]], collapse = "/")
      )

      objfile <- file.path(inFolder, "txi.RDS")

      if (file.exists(objfile)) {

        message("txi.RDS detected ✔")

        rv$app_state <- "ready"
        rv$projFolderFull <- inFolder
        rv$colData <- readRDS(file.path(rv$projFolderFull, "colData.RDS"))
        rv$referenceGenomeChoice <- readRDS(file.path(rv$projFolderFull, "referenceGenomeChoice.RDS"))

        if(rv$referenceGenomeChoice!=1){
          rv$OrgDeeBee <- org.Mm.eg.db::org.Mm.eg.db
        } else {
          rv$OrgDeeBee <- org.Hs.eg.db::org.Hs.eg.db
        }

        withProgress(message = "Opening project...", value = 0, {

          incProgress(0.15, detail = "Reading txi object")
          rv$txi <- readRDS(objfile)


          txi_tpm_file <- file.path(rv$projFolderFull, "txi_tpms.RDS")
          if(file.exists(txi_tpm_file)){
            incProgress(0.2, detail = "Reading TPM matrices")
            rv$txi_tpms <- readRDS(txi_tpm_file)
          } else {
            incProgress(0.2, detail = "TPM matrix not detected, calculating and saving")
            rv$txi_tpms <- calculate_txi_tpm(rv$txi, rv$OrgDeeBee, rv$colData)
            saveRDS(rv$txi_tpms, file = txi_tpm_file)
          }

          txi_deseq_file <- file.path(rv$projFolderFull, "txi_deseq.RDS")
          if(file.exists(txi_deseq_file)){
            incProgress(0.2, detail = "Reading txi_deseq.RDS")
            rv$txi_deseq <- readRDS(txi_deseq_file)
          } else {
            incProgress(0.2, detail = "txi_deseq not detected, calculating and saving")
            rv$txi_deseq <- DESeq2::DESeqDataSetFromTximport(rv$txi,
                                                             colData = rv$colData,
                                                             design = ~Group)
            saveRDS(rv$txi_deseq, file = txi_deseq_file)
          }



          txi_deseq_deseq_file <- file.path(rv$projFolderFull, "txi_deseq_deseq.RDS")
          if(file.exists(txi_deseq_deseq_file)){
            incProgress(0.2, detail = "Reading txi_deseq_deseq.RDS")
            rv$txi_deseq_deseq <- readRDS(txi_deseq_deseq_file)
          } else {
            incProgress(0.2, detail = "txi_deseq_deseq not detected, calculating and saving")
            rv$txi_deseq_deseq <- calculate_txi_deseq_deseq(rv$txi_deseq, rv$colData)
            saveRDS(rv$txi_deseq_deseq, file = txi_deseq_deseq_file)
          }

          res_txi_deseq_file <- file.path(rv$projFolderFull, "res_txi_deseq.RDS")
          if(file.exists(res_txi_deseq_file)){
            incProgress(0.35, detail = "Reading DESeq2 results")
            rv$res_txi_deseq <- readRDS(res_txi_deseq_file)
          } else {
            incProgress(0.35, detail = "DESeq2 not detected, calculating and saving")
            rv$res_txi_deseq <- calculate_res_txi_deseq(rv$txi_deseq_deseq, rv$OrgDeeBee)
              DESeq2::DESeqDataSetFromTximport(rv$txi,
                                                             colData = rv$colData,
                                                             design = ~Group)
            saveRDS(rv$res_txi_deseq, file = res_txi_deseq_file)
          }

          res_DEGs_txi_deseq_file <- file.path(rv$projFolderFull, "res_DEGs_txi_deseq.RDS")
          if(file.exists(res_DEGs_txi_deseq_file)){
            incProgress(0.2, detail = "Reading res_DEGs_txi_deseq")
            rv$res_DEGs_txi_deseq <- readRDS(res_DEGs_txi_deseq_file)
          } else {
            incProgress(0.2, detail = "Calculating res_DEGs_txi_deseq and saving")
            rv$res_DEGs_txi_deseq <- calculate_res_DEGs_txi_deseq(rv$res_txi_deseq)
            saveRDS(rv$res_DEGs_txi_deseq, file = res_DEGs_txi_deseq_file)
          }

          txi_deseq_deseq_file <- file.path(rv$projFolderFull, "txi_deseq_deseq.RDS")
          if(file.exists(txi_deseq_deseq_file)){
            incProgress(0.2, detail = "Reading txi_deseq_deseq")
            rv$txi_deseq_deseq <- readRDS(txi_deseq_deseq_file)
          } else {
            incProgress(0.2, detail = "txi_deseq_deseq not detected, calculating and saving")
            rv$txi_deseq_deseq <- calculate_txi_deseq_deseq_file(rv$txi_deseq, rv$colData)
            saveRDS(rv$txi_deseq_deseq, file = txi_deseq_deseq_file)
          }

          incProgress(0.15, detail = "Reading GO enrichments...")

          GO_result_file <- file.path(rv$projFolderFull, "GO_result.RDS")
          if(file.exists(GO_result_file)){
            incProgress(0.2, detail = "Reading GO enrichment results data")
            rv$GO_result <- readRDS(GO_result_file)
          } else {
            incProgress(0.2, detail = "GO not detected; Calculating GO enrichment results and saving")
            rv$GO_result <- calculate_GO_result(rv$res_DEGs_txi_deseq, rv$OrgDeeBee)
            saveRDS(rv$GO_result, file = GO_result_file)
          }

          vst_file <- file.path(rv$projFolderFull, "vst_data.RDS")
          if(file.exists(vst_file)){
            incProgress(0.2, detail = "Reading VST transformed data")
            rv$vst_data <- readRDS(vst_file)
          } else {
            incProgress(0.2, detail = "VST transforming txi data and saving")
            rv$vst_data <- DESeq2::vst(rv$txi_deseq)
            saveRDS(rv$vst_data, file = vst_file)
          }
          incProgress(0.1, detail = "Finished loading")

        })

        if(length(unique(rv$colData$Group))  > 2) {
          rv$multiple_groups <- 1
        } else {
          rv$multiple_groups <- 0
        }
      } else {
        rv$app_state <- "idle"
      }
    })

    output$previousprojname <- renderText({

      req(input$infolder)

      inFolder <- paste0(
        house_path,
        paste0(input$infolder[[1]], collapse = "/")
      )

      objfile <- file.path(inFolder, "txi.RDS")

      if (file.exists(objfile)) {
        "Loaded project detected at the selected destination"
      } else {
        "No project results detected at selected destination."
      }

    })

  })
}
