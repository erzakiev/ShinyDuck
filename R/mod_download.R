mod_download_ui <- function(id) {
  ns <- NS(id)
  tagList(

    tags$script(HTML("
Shiny.addCustomMessageHandler('check_all_boxes', function(value) {

  document.querySelectorAll('.file-select').forEach(function(el){
    el.checked = value;
  });

});
")),
  h3("Download data?"),
  fluidRow(
    column(2,
           actionButton(ns("check_all"), "Check all")
    ),
    column(2,
           actionButton(ns("uncheck_all"), "Uncheck all")
    ),
    column(1,
           checkboxInput(ns("rds_include"), "Include development files (.RDS)?", value = FALSE)
    )
  ),
  DT::DTOutput(ns("download_table")),
  downloadButton(ns("download_selected"), "Zip and Download")
  )
}

mod_download_server <- function(id, rv) {

  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    table_data <- reactive({

      req(rv$projFolderFull)

      include_rds <- isTRUE(input$rds_include)


      req(rv$app_state == "ready")

      if (input$rds_include) {
        filelist <- list.files(rv$projFolderFull, pattern="\\.RDS$|\\.xlsx$|\\.png$")
      } else {
        filelist <- list.files(rv$projFolderFull, pattern="\\.xlsx$|\\.png$")
      }

      sizes <- file.info(file.path(rv$projFolderFull, filelist))$size
      sizes <- round(sizes / 1024^2, 2)

      data.frame(
        selected = '<input type="checkbox" class="file-select" checked>',
        file = filelist,
        size_MB = sizes,
        description = "",
        stringsAsFactors = FALSE
      )

    })



    output$download_table <- DT::renderDT({
      print("renderDT triggered")
      print(head(table_data()))

      DT::datatable(
        table_data(),
        escape = FALSE,
        selection = "none",
        rownames = FALSE,
        options = list(
          dom = 't',
          pageLength = 50
        ),
        callback = htmlwidgets::JS(
          paste0(
            "
table.on('change', '.file-select', function() {

  var checked = [];

  table.$('.file-select').each(function(i) {

    if(this.checked){
      checked.push(i+1);
    }

  });

  Shiny.setInputValue('", ns("checked_rows"), "', checked, {priority:'event'});

});
"
          )
        )
      )
    })


    observeEvent(input$check_all, {

      session$sendCustomMessage("check_all_boxes", TRUE)

    })

    observeEvent(input$uncheck_all, {

      session$sendCustomMessage("check_all_boxes", FALSE)

    })


    output$download_selected <- downloadHandler(

      filename = function() {
        paste0("ShinyDucks_results_", basename(rv$projFolderFull), ".zip")
      },

      content = function(file) {

        rows <- input$checked_rows
        if (is.null(rows)) rows <- seq_len(nrow(table_data()))

        selected <- table_data()$file[rows]

        withProgress(message = "Preparing archive...", value = 0, {

          incProgress(0.3, detail = paste("Archiving", length(selected), "files"))

          zip::zip(
            zipfile = file,
            files = selected,
            root = rv$projFolderFull
          )

          incProgress(0.7)


        })

      }
    )
    })
}


