mod_tpm_table_ui <- function(id) {
  ns <- NS(id)
tagList(
  h3("Transcripts-per-million (TPM) values:"),
  #fluidRow(
  #  column(
  #    3,
  #
  #    shinyWidgets::prettySwitch(
  #      inputId = ns("color_scale_mode"),
  #      label = "Global color scaling",
  #      value = TRUE
  #    )),
  #  column(
  #    3,
  #
  #    shinyWidgets::prettySwitch(
  #      inputId = ns("log_transform"),
  #      label = "log10(TPM+1) (= like in real heatmaps)",
  #      value = FALSE
  #    ))),
  DT::DTOutput(ns("TPM"))
)
}

mod_tpm_table_server <- function(id, rv) {

  moduleServer(id, function(input, output, session) {

    output$TPM <- DT::renderDT({

      #df <- as.data.frame(rv$txi_tpms)
      #df_subset <- t(apply(df[2:nrow(df),3:ncol(df)], MARGIN = 1, as.numeric))
      #df[2:nrow(df),3:ncol(df)] <- format_df_numbers(as.data.frame(df_subset), digits = 3)


      #df <- df[-1,]
      #df_subset <- t(apply(df[,3:ncol(df)], MARGIN = 1, as.numeric))
      #df[,3:ncol(df)] <- format_df_numbers(as.data.frame(df_subset), digits = 3)

      df <- as.data.frame(rv$txi$abundance)
      df <- convertDfRownamesEns2Symb(df, rv$OrgDeeBee)

      ord <- order(rv$colData$Group)
      df <- df[,ord]

      #colnames(df) <- gsub('.fastq.trimmed', replacement = '', colnames(df))
      #colnames(df) <- gsub('_001', replacement = '', colnames(df))

      df <- cbind(SYMBOL=rownames(df), df)

      container <- make_TPM_container(
        df,
        rv$colData[ord,]
      )

      DT::datatable(
        df,
        extensions = "Buttons",
        options = list(
          dom = "Blfrtip",
          buttons =
            list(list(extend='copy'),
                 list(extend='csv', filename = paste0("TPM_table")),
                 list(extend='excel', filename = paste0("TPM_table")),
                 list(extend='pdf', filename = paste0("TPM_table")),
                 list(extend='print')),
          pageLength = 25,
          scrollX = TRUE
        ),
        container = container,
        rownames = F
      ) |> color_tpm_table(
        df,
        #global_scale = input$color_scale_mode,
        #log_scale = input$log_transform
        global_scale = T,
        log_scale = F
      ) |>
        add_group_separators(
          df,
          rv$colData
        )
    }, server = FALSE)

  })
}
