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
      "Run analysis",
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

    observeEvent(input$referenceGenomeChoice, {
      rv$referenceGenomeChoice <- input$referenceGenomeChoice
    })

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
        rv$txi <- readRDS(objfile)
        rv$txi_tpms <- readRDS(file.path(inFolder, "txi_tpms.RDS"))
        rv$txi_deseq <- readRDS(file.path(inFolder, "txi_deseq.RDS"))
        rv$res_txi_deseq <- readRDS(file.path(inFolder, "res_txi_deseq.RDS"))
        rv$res_DEGs_txi_deseq <- readRDS(file.path(inFolder, "res_DEGs_txi_deseq.RDS"))
        rv$txi_deseq_deseq <- readRDS(file.path(inFolder, "txi_deseq_deseq.RDS"))
        rv$colData <- readRDS(file.path(inFolder, "colData.RDS"))
          if(rv$colData$Group %>% nlevels() > 2) {
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
        "Loading project detected at selected destination."
      } else {
        "No project results detected at selected destination."
      }

    })

  })
}
