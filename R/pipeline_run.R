library(shiny)
library(shinyFiles)
library(DT)

# -------------------------------------------------------------------------
# UI
# -------------------------------------------------------------------------

mod_rnaseq_ui <- function(id) {
  
  ns <- NS(id)
  
  tagList(
    
    h3("nf-core/rnaseq"),
    
    shinyFilesButton(
      id = ns("fastq"),
      label = "Select FASTQ files",
      title = "Select FASTQ files",
      multiple = TRUE
    ),
    
    br(),
    br(),
    
    h4("Selected FASTQ files"),
    
    DTOutput(ns("fastq_table")),
    
    br(),
    
    h4("Samplesheet"),
    
    DTOutput(ns("samplesheet")),
    
    br(),
    
    downloadButton(
      ns("download_samplesheet"),
      "Download samplesheet.csv"
    )
    
  )
}


# -------------------------------------------------------------------------
# Server
# -------------------------------------------------------------------------

mod_rnaseq_server <- function(
    id,
    fastq_root
) {
  
  moduleServer(id, function(input, output, session) {
    
    ns <- session$ns
    
    # ---------------------------------------------------------------
    # FASTQ file selector
    # ---------------------------------------------------------------
    
    roots <- c(
      FASTQ = fastq_root
    )
    
    shinyFileChoose(
      input,
      "fastq",
      roots = roots,
      session = session,
      filetypes = c(
        "fastq",
        "fq",
        "gz"
      )
    )
    
    
    # ---------------------------------------------------------------
    # Selected files
    # ---------------------------------------------------------------
    
    selected_files <- reactive({
      
      req(input$fastq)
      
      files <- parseFilePaths(
        roots,
        input$fastq
      )
      
      files$datapath
      
    })
    
    
    # ---------------------------------------------------------------
    # Display selected files
    # ---------------------------------------------------------------
    
    output$fastq_table <- renderDT({
      
      files <- selected_files()
      
      datatable(
        data.frame(
          file = basename(files),
          path = files
        ),
        rownames = FALSE,
        options = list(
          pageLength = 10
        )
      )
      
    })
    
    
    # ---------------------------------------------------------------
    # Detect R1 / R2
    # ---------------------------------------------------------------
    
    samplesheet <- reactive({
      
      files <- selected_files()
      
      filenames <- basename(files)
      
      is_r1 <- grepl(
        "(_R?1|_R1_001)\\.f(ast)?q\\.gz$",
        filenames,
        ignore.case = TRUE
      )
      
      is_r2 <- grepl(
        "(_R?2|_R2_001)\\.f(ast)?q\\.gz$",
        filenames,
        ignore.case = TRUE
      )
      
      r1 <- files[is_r1]
      r2 <- files[is_r2]
      
      if (length(r1) == 0) {
        return(data.frame())
      }
      
      # Remove R1 suffix to obtain sample name
      get_sample_name <- function(x) {
        
        x <- basename(x)
        
        x <- sub(
          "_R1(_001)?\\.f(ast)?q\\.gz$",
          "",
          x,
          ignore.case = TRUE
        )
        
        x <- sub(
          "_R2(_001)?\\.f(ast)?q\\.gz$",
          "",
          x,
          ignore.case = TRUE
        )
        
        x
      }
      
      samples_r1 <- get_sample_name(r1)
      samples_r2 <- get_sample_name(r2)
      
      sample_names <- unique(
        c(samples_r1, samples_r2)
      )
      
      result <- data.frame(
        sample = sample_names,
        fastq_1 = NA_character_,
        fastq_2 = NA_character_,
        strandedness = "auto",
        stringsAsFactors = FALSE
      )
      
      for (i in seq_along(sample_names)) {
        
        s <- sample_names[i]
        
        idx1 <- which(samples_r1 == s)
        
        if (length(idx1) > 0) {
          result$fastq_1[i] <- r1[idx1[1]]
        }
        
        idx2 <- which(samples_r2 == s)
        
        if (length(idx2) > 0) {
          result$fastq_2[i] <- r2[idx2[1]]
        }
      }
      
      result
      
    })
    
    
    # ---------------------------------------------------------------
    # Display samplesheet
    # ---------------------------------------------------------------
    
    output$samplesheet <- renderDT({
      
      datatable(
        samplesheet(),
        editable = TRUE,
        rownames = FALSE,
        options = list(
          pageLength = 10,
          scrollX = TRUE
        )
      )
      
    })
    
    
    # ---------------------------------------------------------------
    # Download samplesheet
    # ---------------------------------------------------------------
    
    output$download_samplesheet <- downloadHandler(
      
      filename = function() {
        "samplesheet.csv"
      },
      
      content = function(file) {
        
        write.csv(
          samplesheet(),
          file,
          row.names = FALSE,
          quote = FALSE
        )
        
      }
      
    )
    
  })
}