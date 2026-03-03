mod_tpm_table_ui <- function(id) {
  ns <- NS(id)
tagList(
  textOutput(ns("text_TPM")),
  DT::DTOutput(ns("TPM"))
)
}

mod_tpm_table_server <- function(id, rv) {

  moduleServer(id, function(input, output, session) {

    output$text_TPM <- renderText({
      "Transcripts-per-million (TPM) values:"
    })

    output$TPM <- DT::renderDT({

      DT::datatable(
        as.data.frame(rv$txi_tpms),
        extensions = "Buttons",
        options = list(
          dom = "Blfrtip",
          buttons = c("copy", "csv", "excel", "pdf", "print"),
          pageLength = 25,
          scrollX = TRUE
        ),
        rownames = FALSE
      )

    }, server = FALSE)

  })
}
