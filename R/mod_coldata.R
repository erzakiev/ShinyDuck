mod_coldata_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(
        width = 6,
        style = "border-right: 5px solid #3d8cbc; padding-right: 20px;",
        DT::DTOutput(ns("GroupsPrompt")),
        actionButton(ns("groups_rebase"),
                     label = "Modify background factor / groups?",
                     icon = shiny::icon("rotate"))
      ),
      column(
        width = 6,
        h3('Subset into another folder?'),
        textOutput(ns("GroupsPromptTextFork")),
        textInput(ns("folder2fork2"),
                  label = "Please specify new folder name"),
        uiOutput(ns("ForkingPrompt")),
        uiOutput(ns("ForkingPrompt2")),
        actionButton(ns("validateforking"),
                     label = "Validate choice"),
        verbatimTextOutput(ns("forkingfeedback"))
      )
    )
  )
}


mod_coldata_server <- function(id, rv) {

  moduleServer(id, function(input, output, session) {


    ns <- session$ns

    # ---- TABLE RENDER ----
    output$GroupsPrompt <- DT::renderDT({

      message("Rendering GroupsPrompt")

      df <- rv$colData

      df$Control <- sprintf(
        '<input type="radio" name="control_radio" value="%s" %s/>',
        1:nrow(df),
        ifelse(df$Control, "checked", "")
      )

      DT::datatable(
        df,
        escape = FALSE,
        editable = list(
          target = "cell",
          disable = list(columns = c(0,2))
        ),
        rownames = FALSE,
        options = list(
          pageLength = 20,
          dom = "tip"
        ),

        callback = htmlwidgets::JS(sprintf("
  table.on('change', 'input[type=\"radio\"]', function() {
    Shiny.setInputValue('%s', this.value, {priority: 'event'});
  });
", session$ns("control_selected")))

      )

    })

    observeEvent(input$control_selected, {

      i <- as.numeric(input$control_selected)
      message(paste0('selected ', i, ' of as.numeric(input$control_selected)'))

      rv$colData$Control <- FALSE
      rv$colData$Control[i] <- TRUE
      saveRDS(rv$colData, file = file.path(rv$projFolderFull, "colData.RDS"))
    })

    # ---- CELL EDIT HANDLER ----
    observeEvent(input$GroupsPrompt_cell_edit, {

      info <- input$GroupsPrompt_cell_edit

      i <- info$row
      j <- info$col + 1 # important lol indexing in DT is different from the R's one
      v <- info$value

      rv$colData[i, j] <- v

      # update multi-group flag
      if ("Group" %in% colnames(rv$colData)) {
        rv$multiple_groups <- ifelse(length(unique(rv$colData$Group)) > 2, yes = 1, no = 0)
      }
      rv$colData[,'Control'] <- as.logical(rv$colData[,'Control'])

    })

    # ---- SAVE GROUPS ----
    observeEvent(input$groups_specified, {

      req(rv$colData)
      req(rv$projFolderFull)

      write.table(
        rv$colData,
        file = file.path(rv$projFolderFull, "colData.tsv"),
        sep = "\t",
        row.names = FALSE
      )

      saveRDS(
        rv$colData,
        file = file.path(rv$projFolderFull, "colData.RDS")
      )

      message("colData saved")

    })

    # ---- REBASE ----

    observeEvent(input$groups_rebase, {
      req(rv$colData)
      req(rv$projFolderFull)


      withProgress(message = "Rebasing...", value = 0, {


        saveRDS(
          rv$colData,
          file = file.path(rv$projFolderFull, "colData.RDS")
        )

        # recalculating all the stuff as the basis of that
        reff <- as.character(rv$colData[which(rv$colData[,3]),2][1])
        rv$colData$Group <- as.factor(rv$colData$Group)
        if(length(which(rv$colData[,3]))>0) rv$colData$Group <- relevel(x = rv$colData$Group, ref = reff)

        incProgress(0.05, detail = 'Remaking the DESeq2::DESeqDataSet')
        rv$txi_deseq <- DESeq2::DESeqDataSetFromTximport(rv$txi, colData = rv$colData, design = ~Group)

        incProgress(0.05, detail = 'Saving the DESeq2::DESeqDataSet')
        saveRDS(
          rv$txi_deseq,
          file = file.path(rv$projFolderFull, "txi_deseq.RDS")
        )

        incProgress(0.05, detail = 'Filtering lowly expressed genes before DEGs')


        incProgress(0.05, detail = 'Recalculating txi_deseq_deseq')
        rv$txi_deseq_deseq <- calculate_txi_deseq_deseq(rv$txi_deseq, rv$colData)

        incProgress(0.05, detail = 'Saving txi_deseq_deseq')
        saveRDS(rv$txi_deseq_deseq,
                file = file.path(rv$projFolderFull, 'txi_deseq_deseq.RDS'))


        rv$res_txi_deseq <- calculate_res_txi_deseq(rv$txi_deseq_deseq, rv$OrgDeeBee)
        saveRDS(rv$res_txi_deseq,
                file = file.path(rv$projFolderFull, 'res_txi_deseq.RDS'))
        #openxlsx::write.xlsx(rv$res_txi_deseq, file = file.path(rv$projFolderFull,'DEGs_full.xlsx'))

        rv$res_DEGs_txi_deseq <- calculate_res_DEGs_txi_deseq(rv$res_txi_deseq)
        saveRDS(rv$res_DEGs_txi_deseq,
                file = file.path(rv$projFolderFull, 'res_DEGs_txi_deseq.RDS'))

        #if(class(rv$res_DEGs_txi_deseq)=="data.frame"){
        #  openxlsx::write.xlsx(as.data.frame(rv$res_DEGs_txi_deseq), file = file.path(rv$projFolderFull,'DEGs.xlsx'))
        #} else {
        #  openxlsx::write.xlsx(rv$res_DEGs_txi_deseq, file = file.path(rv$projFolderFull,'DEGs.xlsx'))
        #}

        incProgress(0.05, detail = 'Recalculating and saving GO_result')
        rv$GO_result <- calculate_GO_result(rv$res_DEGs_txi_deseq, rv$OrgDeeBee)
        saveRDS(
          rv$GO_result,
          file = file.path(rv$projFolderFull, "GO_result.RDS")
        )
        #openxlsx::write.xlsx(toWrite, file = file.path(rv$projFolderFull,'GOs.xlsx'))

        incProgress(0.05, detail = 'Recalculating and saving vst_data')
        rv$vst_data <- DESeq2::vst(rv$txi_deseq)
        saveRDS(
          rv$vst_data,
          file = file.path(rv$projFolderFull, "vst_data.RDS")
        )
        message("colData rebased")


      })
    })

    # ---- FORKING UI ----
    output$ForkingPrompt <- renderUI({

      req(rv$app_state == "ready")
      req(rv$colData)

      shiny::selectInput(
        ns("samples2fork"),
        label = "Samples to fork into another directory",
        choices = rv$colData$Sample,
        multiple = TRUE
      )
    })


    output$ForkingPrompt2 <- renderUI({

      req(rv$app_state == "ready")
      req(rv$colData)

      shiny::selectInput(
        ns("groups2fork"),
        label = "... or select whole groups to fork into another directory",
        choices = unique(rv$colData$Group),
        multiple = TRUE
      )

    })

    # ---- FORKING LOGIC ----
    observeEvent(input$validateforking, {

      req(input$folder2fork2)
      req(input$samples2fork)

      foldr <- gsub("[^A-Za-z0-9_.]", "_", input$folder2fork2)

      subst <- which(rv$colData$Sample %in% input$samples2fork)

      colData_subset <- rv$colData[subst, , drop = FALSE]

      txi_subset <- rv$txi
      txi_subset$abundance <- txi_subset$abundance[, subst]
      txi_subset$counts    <- txi_subset$counts[, subst]
      txi_subset$length    <- txi_subset$length[, subst]

      colData_subset$Control <- FALSE
      colData_subset$Control[1] <- TRUE

      new_dir <- file.path(dirname(rv$projFolderFull), foldr)
      dir.create(new_dir, showWarnings = FALSE)

      withProgress(message = "Forking...", value = 0, {
        incProgress(0.05, detail = 'saving new txi object')
        saveRDS(txi_subset, file.path(new_dir, "txi.RDS"))
        saveRDS(colData_subset, file.path(new_dir, "colData.RDS"))
        saveRDS(rv$referenceGenomeChoice,
                file.path(new_dir, "referenceGenomeChoice.RDS"))

        reff <- as.character(colData_subset[which(colData_subset[,3]),2][1])
        colData_subset$Group <- as.factor(colData_subset$Group)
        if(length(which(colData_subset[,3]))>0) colData_subset$Group <- relevel(x = colData_subset$Group, ref = reff)
        incProgress(0.05, detail = 'recreating and saving the new txi_deseq object')
        txi_deseq <- DESeq2::DESeqDataSetFromTximport(txi_subset, colData = colData_subset, design = ~Group)

        saveRDS(
          txi_deseq,
          file = file.path(new_dir, "txi_deseq.RDS")
        )

        saveRDS(
          calculate_txi_tpm(txi = txi_subset,
                            #colData = colData_subset,
                            orgdb = rv$OrgDeeBee),
          file = file.path(new_dir, "txi_tpms.RDS")
        )

        incProgress(0.05, detail = 'calculating and saving the new txi_deseq_deseq object')
        txi_deseq_deseq <- calculate_txi_deseq_deseq(txi_deseq, colData_subset)
        saveRDS(txi_deseq_deseq,
                file = file.path(new_dir, 'txi_deseq_deseq.RDS'))

        incProgress(0.05, detail = 'calculating and saving the new res_txi_deseq object')
        res_txi_deseq <- calculate_res_txi_deseq(txi_deseq_deseq, rv$OrgDeeBee)
        saveRDS(res_txi_deseq,
                file = file.path(new_dir, 'res_txi_deseq.RDS'))
        #openxlsx::write.xlsx(res_txi_deseq, file = file.path(new_dir,'DEGs_full.xlsx'))

        incProgress(0.05, detail = 'calculating and saving the new res_DEGs_txi_deseq object')
        res_DEGs_txi_deseq <- calculate_res_DEGs_txi_deseq(res_txi_deseq)
        saveRDS(res_DEGs_txi_deseq,
                file = file.path(new_dir, 'res_DEGs_txi_deseq.RDS'))

        #if(class(toRet)=="data.frame"){
       #   openxlsx::write.xlsx(as.data.frame(toRet), file = file.path(new_dir,'DEGs.xlsx'))
        #} else {
        #  openxlsx::write.xlsx(toRet, file = file.path(new_dir,'DEGs.xlsx'))
        #}

        # recalculating GO_result

        incProgress(0.05, detail = 'calculating and saving the new GO_result object')
        GO_result <- calculate_GO_result(res_DEGs_txi_deseq,
                                         rv$OrgDeeBee)
        saveRDS(GO_result,
                file = file.path(new_dir, 'GO_result.RDS'))
        #openxlsx::write.xlsx(reshape_GO_result_for_xlsx(GO_result), file = file.path(new_dir,'GOs.xlsx'))

        incProgress(0.5, detail = 'calculating and saving the new vst_data object')
        saveRDS(
          DESeq2::vst(txi_deseq),
          file = file.path(new_dir, "vst_data.RDS")
        )
      })

      output$forkingfeedback <- renderText(
        paste0("Successfully copied samples into folder ", foldr)
      )

    })


    observeEvent(input$validateforking, {

      req(input$folder2fork2)
      req(input$groups2fork)

      foldr <- gsub("[^A-Za-z0-9_.]", "_", input$folder2fork2)

      subst <- which(rv$colData$Group %in% input$groups2fork)

      colData_subset <- rv$colData[subst, , drop = FALSE]

      txi_subset <- rv$txi
      txi_subset$abundance <- txi_subset$abundance[, subst]
      txi_subset$counts    <- txi_subset$counts[, subst]
      txi_subset$length    <- txi_subset$length[, subst]

      new_dir <- file.path(dirname(rv$projFolderFull), foldr)
      dir.create(new_dir, showWarnings = FALSE)
      colData_subset$Control <- FALSE
      colData_subset$Control[1] <- TRUE

      withProgress(message = "Forking...", value = 0, {
        incProgress(0.05, detail = 'saving new txi object')
        saveRDS(txi_subset, file.path(new_dir, "txi.RDS"))
        saveRDS(colData_subset, file.path(new_dir, "colData.RDS"))
        saveRDS(rv$referenceGenomeChoice,
              file.path(new_dir, "referenceGenomeChoice.RDS"))

      reff <- as.character(colData_subset[which(colData_subset[,3]),2][1])
      colData_subset$Group <- as.factor(colData_subset$Group)
      if(length(which(colData_subset[,3]))>0) colData_subset$Group <- relevel(x = colData_subset$Group, ref = reff)
      incProgress(0.05, detail = 'recreating and saving the new txi_deseq object')
      txi_deseq <- DESeq2::DESeqDataSetFromTximport(txi_subset, colData = colData_subset, design = ~Group)

      saveRDS(
        txi_deseq,
        file = file.path(new_dir, "txi_deseq.RDS")
      )

    saveRDS(
      calculate_txi_tpm(txi = txi_subset,
                        #colData = colData_subset,
                        orgdb = rv$OrgDeeBee),
      file = file.path(new_dir, "txi_tpms.RDS")
    )

    incProgress(0.05, detail = 'calculating and saving the new txi_deseq_deseq object')
     txi_deseq_deseq <- calculate_txi_deseq_deseq(txi_deseq, colData_subset)
     saveRDS(txi_deseq_deseq,
              file = file.path(new_dir, 'txi_deseq_deseq.RDS'))

     incProgress(0.05, detail = 'calculating and saving the new res_txi_deseq object')
      res_txi_deseq <- calculate_res_txi_deseq(txi_deseq_deseq, rv$OrgDeeBee)
      saveRDS(res_txi_deseq,
              file = file.path(new_dir, 'res_txi_deseq.RDS'))
      #openxlsx::write.xlsx(res_txi_deseq, file = file.path(new_dir,'DEGs_full.xlsx'))

      incProgress(0.05, detail = 'calculating and saving the new res_DEGs_txi_deseq object')
      res_DEGs_txi_deseq <- calculate_res_DEGs_txi_deseq(res_txi_deseq)
      saveRDS(res_DEGs_txi_deseq,
              file = file.path(new_dir, 'res_DEGs_txi_deseq.RDS'))

      #if(class(toRet)=="data.frame"){
        #openxlsx::write.xlsx(as.data.frame(toRet), file = file.path(new_dir,'DEGs.xlsx'))
      #} else {
        #openxlsx::write.xlsx(toRet, file = file.path(new_dir,'DEGs.xlsx'))
      #}

      # recalculating GO_result

      incProgress(0.05, detail = 'calculating and saving the new GO_result object')
      GO_result <- calculate_GO_result(res_DEGs_txi_deseq,
                                       rv$OrgDeeBee)
      saveRDS(GO_result,
              file = file.path(new_dir, 'GO_result.RDS'))
      #openxlsx::write.xlsx(reshape_GO_result_for_xlsx(GO_result), file = file.path(new_dir,'GOs.xlsx'))

      incProgress(0.5, detail = 'calculating and saving the new vst_data object')
      saveRDS(
        DESeq2::vst(txi_deseq),
        file = file.path(new_dir, "vst_data.RDS")
      )
    })


      output$forkingfeedback <- renderText(
        paste0("Successfully copied samples into folder ", foldr)
      )

    })

  })
}
