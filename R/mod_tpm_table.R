mod_tpm_table_ui <- function(id) {
  ns <- NS(id)
tagList(
  h3("Transcripts-per-million (TPM) values:"),
  DT::DTOutput(ns("TPM"))
)
}

mod_tpm_table_server <- function(id, rv) {

  moduleServer(id, function(input, output, session) {

    output$TPM <- DT::renderDT({

      df <- as.data.frame(rv$txi_tpms)
      df_subset <- t(apply(df[2:nrow(df),3:ncol(df)], MARGIN = 1, as.numeric))

      df[2:nrow(df),3:ncol(df)] <- format_df_numbers(as.data.frame(df_subset), digits = 3)
      colnames(df) <- gsub('.fastq.trimmed', replacement = '', colnames(df))
      colnames(df) <- gsub('_001', replacement = '', colnames(df))

      DT::datatable(
        df,
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
