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

      saveRDS(
        rv$colData,
        file = file.path(rv$projFolderFull, "colData.RDS")
      )

      # recalculating all the stuff as the basis of that
      reff <- as.character(rv$colData[which(rv$colData[,3]),2][1])
      if(length(which(rv$colData[,3]))>0) rv$colData$Group <- relevel(x = rv$colData$Group, ref = reff)
      rv$txi_deseq <- DESeq2::DESeqDataSetFromTximport(rv$txi, colData = rv$colData, design = ~Group)

      saveRDS(
        rv$txi_deseq,
        file = file.path(rv$projFolderFull, "txi_deseq.RDS")
      )

      matr <- rv$txi_deseq
      smallestGroupSize <- floor(min(table(rv$colData$Group)))
      keep <- rowSums(DESeq2::counts(matr) >= 10) >= smallestGroupSize
      matr <- matr[keep,]
      toRet <- DESeq2::DESeq(matr)
      rv$txi_deseq_deseq <- toRet

      saveRDS(rv$txi_deseq_deseq,
              file = file.path(rv$projFolderFull, 'txi_deseq_deseq.RDS'))

      toRet <- toRet2 <- toRet3 <- toWrite <- tytl <- list()
      if(length(DESeq2::resultsNames(rv$txi_deseq_deseq))<=2){
        toRet <- DESeq2::results(rv$txi_deseq_deseq)
        toRet$log10padj <- log10(toRet$padj)
        toRet$Significance <- 'Non-significant'
        toRet$Significance[toRet$padj<0.05] <- 'Significant (padj < 0.05)'

        toRet2 <- data.frame(ENSEMBL=rownames(toRet),
                             toRet[,1:2],
                             abs_log2FoldChange=abs(toRet[,2]),toRet[,3:ncol(toRet)])
        annots <- AnnotationDbi::select(rv$OrgDeeBee, keys=rownames(toRet2),
                                        columns="SYMBOL", keytype="ENSEMBL")

        toRet3 <- merge(annots, toRet2, by.x="ENSEMBL", by.y="ENSEMBL")
        toRet3 <- toRet3[order(toRet3$padj),]

      }  else {
        toRet3 <- toRet <- list()
        for (i in 2:length(DESeq2::resultsNames(rv$txi_deseq_deseq))){
          toRet[[(i-1)]] <- DESeq2::results(rv$txi_deseq_deseq, name = DESeq2::resultsNames(rv$txi_deseq_deseq)[i] )
          toRet[[(i-1)]]$log10padj <- log10( toRet[[(i-1)]]$padj)
          toRet[[(i-1)]]$Significance <- 'Non-significant'
          toRet[[(i-1)]]$Significance[ toRet[[(i-1)]]$padj<0.05] <- 'Significant (padj < 0.05)'

          toRet2 <- data.frame(ENSEMBL=rownames(toRet[[(i-1)]]),
                               toRet[[(i-1)]][,1:2],
                               abs_log2FoldChange=abs(toRet[[(i-1)]][,2]),toRet[[(i-1)]][,3:ncol(toRet[[(i-1)]])])

          annots <- AnnotationDbi::select(rv$OrgDeeBee, keys=rownames(toRet2),
                                          columns="SYMBOL", keytype="ENSEMBL")

          toRet3[[(i-1)]] <- merge(annots,
                                   toRet2,
                                   by.x="ENSEMBL",
                                   by.y="ENSEMBL")
          toRet3[[(i-1)]] <- toRet3[[(i-1)]][order(toRet3[[(i-1)]]$padj),]

        }
        names(toRet3) <- DESeq2::resultsNames(rv$txi_deseq_deseq)[2:length(DESeq2::resultsNames(rv$txi_deseq_deseq))]
      }
      rv$res_txi_deseq <- toRet3
      saveRDS(
        rv$res_txi_deseq,
        file = file.path(rv$projFolderFull,'res_txi_deseq.RDS')
      )
      openxlsx::write.xlsx(rv$res_txi_deseq,
                           file = file.path(rv$projFolderFull,'DEGs_full.xlsx'))

      # recalculating res_DEGs_txi_deseq
      toRet <- rv$res_txi_deseq
      if(class(toRet)=="data.frame"){
        toRet <- toRet[which(toRet$padj<0.05),]
        if(nrow(toRet)==0) return({as.data.frame("oops, none significant")})
        toRet <- toRet[order(toRet$padj),]
        openxlsx::write.xlsx(as.data.frame(toRet), file = file.path(rv$projFolderFull,'DEGs.xlsx'))
      } else {
        for (i in 1:length(toRet)){
          toRet[[i]] <- toRet[[i]][which(toRet[[i]]$padj<0.05),]
          if(nrow(toRet[[i]])==0){
            toRet[[i]] <- "oops, none significant"
          }
          else {
            toRet[[i]] <- toRet[[i]][order(toRet[[i]]$padj),]
          }
        }
        openxlsx::write.xlsx(toRet, file = file.path(rv$projFolderFull,'DEGs.xlsx'))
      }
      rv$res_DEGs_txi_deseq <- toRet
      saveRDS(
        rv$res_DEGs_txi_deseq,
        file = file.path(rv$projFolderFull,'res_DEGs_txi_deseq.RDS')
      )

      # recalculating GO_result
      toRet <- toRet2 <- toRet3 <- toWrite <- tytl <- list()
      if(class(rv$res_DEGs_txi_deseq)=='list'){
        # multiple groups
        for (i in 1:length(rv$res_DEGs_txi_deseq)){
          tytl[[i]] <- paste(strsplit(names(rv$res_txi_deseq[i]), split = "_")[[1]][2:4], collapse=' ')
          if(is.null(dim(rv$res_DEGs_txi_deseq[[i]]))){
            toRet[[i]] <- NULL
            toWrite[[i]] <- NULL
          } else {
            toRet[[i]] <- clusterProfiler::enrichGO(gene = rv$res_DEGs_txi_deseq[[i]][,2],
                                   keyType = "SYMBOL",
                                   OrgDb = rv$OrgDeeBee,
                                   ont = "BP",
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.05,
                                   readable = TRUE)

            toWrite[[i]] <- as.data.frame(toRet[[i]])
            toWrite[[i]]$geneID <- gsub(pattern='/', replacement=' ', toWrite[[i]]$geneID)
          }
        }
        names(toRet) <- unlist(tytl)
        saveRDS(toRet, file = file.path(rv$projFolderFull,'GO_result.RDS'))
        rv$GO_result <- toRet
        names(toRet) <- gsub(pattern = 'Group ', replacement = '', x = names(toRet))
        names(toRet) <- substr(names(toRet), start = 1, stop = 30)
        openxlsx::write.xlsx(toWrite, file = file.path(rv$projFolderFull,'GOs.xlsx'))
        return(toRet)

      } else {
        # 1 vs 1
        toRet <- clusterProfiler::enrichGO(gene = rv$res_DEGs_txi_deseq[,2],
                          keyType = "SYMBOL",
                          OrgDb = rv$OrgDeeBee,
                          ont = "BP",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05,
                          readable = TRUE)
        saveRDS(toRet, file = file.path(rv$projFolderFull,'GO_result.RDS'))
        rv$GO_result <- toRet
        toWrite <- as.data.frame(toRet)
        toWrite$geneID <- gsub(pattern='/', replacement=' ', toWrite$geneID)
        openxlsx::write.xlsx(toWrite, file = file.path(rv$projFolderFull,'GOs.xlsx'))
        return(toRet)
      }

      saveRDS(
        DESeq2::vst(rv$txi_deseq),
        file = file.path(rv$projFolderFull, "vst_data.RDS")
      )
      message("colData rebased")
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

      saveRDS(txi_subset, file.path(new_dir, "txi.RDS"))
      saveRDS(colData_subset, file.path(new_dir, "colData.RDS"))
      saveRDS(rv$referenceGenomeChoice,
              file.path(new_dir, "referenceGenomeChoice.RDS"))


      reff <- as.character(colData_subset[which(colData_subset[,3]),2][1])
      if(length(which(colData_subset[,3]))>0) colData_subset$Group <- relevel(x = colData_subset$Group, ref = reff)
      txi_deseq <- DESeq2::DESeqDataSetFromTximport(txi_subset, colData = colData_subset, design = ~Group)

      #res <- format_df_numbers(rv$txi$abundance, digits = 3)
      #res <- data.frame(ENSEMBL=rownames(res), res)
      #annots <- AnnotationDbi::select(rv$OrgDeeBee, keys=rownames(res),
      #                                columns="SYMBOL", keytype="ENSEMBL")
      #result <- merge(annots, res, by.x="ENSEMBL", by.y="ENSEMBL")
      #result <- rbind(c('-','-',colData_subset$Group),result)

      #openxlsx::write.xlsx(result, file = file.path(new_dir, 'TPMs.xlsx'))
      #saveRDS(
      #  result,
      #  file = file.path(new_dir, 'txi_tpms.RDS')
      #)

      saveRDS(
        txi_deseq,
        file = file.path(new_dir, "txi_deseq.RDS")
      )

      matr <- txi_deseq
      smallestGroupSize <- floor(min(table(colData_subset$Group)))
      keep <- rowSums(DESeq2::counts(matr) >= 10) >= smallestGroupSize
      matr <- matr[keep,]
      toRet <- DESeq2::DESeq(matr)
      txi_deseq_deseq <- toRet

      saveRDS(txi_deseq_deseq,
              file = file.path(new_dir, 'txi_deseq_deseq.RDS'))

      toRet <- toRet2 <- toRet3 <- toWrite <- tytl <- list()
      if(length(DESeq2::resultsNames(txi_deseq_deseq))<=2){
        toRet <- DESeq2::results(txi_deseq_deseq)
        toRet$log10padj <- log10(toRet$padj)
        toRet$Significance <- 'Non-significant'
        toRet$Significance[toRet$padj<0.05] <- 'Significant (padj < 0.05)'

        toRet2 <- data.frame(ENSEMBL=rownames(toRet),
                             toRet[,1:2],
                             abs_log2FoldChange=abs(toRet[,2]),toRet[,3:ncol(toRet)])
        annots <- AnnotationDbi::select(rv$OrgDeeBee, keys=rownames(toRet2),
                                        columns="SYMBOL", keytype="ENSEMBL")

        toRet3 <- merge(annots, toRet2, by.x="ENSEMBL", by.y="ENSEMBL")
        toRet3 <- toRet3[order(toRet3$padj),]

      }  else {
        toRet3 <- toRet <- list()
        for (i in 2:length(DESeq2::resultsNames(txi_deseq_deseq))){
          toRet[[(i-1)]] <- DESeq2::results(txi_deseq_deseq, name = DESeq2::resultsNames(txi_deseq_deseq)[i] )
          toRet[[(i-1)]]$log10padj <- log10( toRet[[(i-1)]]$padj)
          toRet[[(i-1)]]$Significance <- 'Non-significant'
          toRet[[(i-1)]]$Significance[ toRet[[(i-1)]]$padj<0.05] <- 'Significant (padj < 0.05)'

          toRet2 <- data.frame(ENSEMBL=rownames(toRet[[(i-1)]]),
                               toRet[[(i-1)]][,1:2],
                               abs_log2FoldChange=abs(toRet[[(i-1)]][,2]),toRet[[(i-1)]][,3:ncol(toRet[[(i-1)]])])

          annots <- AnnotationDbi::select(rv$OrgDeeBee, keys=rownames(toRet2),
                                          columns="SYMBOL", keytype="ENSEMBL")

          toRet3[[(i-1)]] <- merge(annots,
                                   toRet2,
                                   by.x="ENSEMBL",
                                   by.y="ENSEMBL")
          toRet3[[(i-1)]] <- toRet3[[(i-1)]][order(toRet3[[(i-1)]]$padj),]

        }
        names(toRet3) <- DESeq2::resultsNames(txi_deseq_deseq)[2:length(DESeq2::resultsNames(txi_deseq_deseq))]
      }
      res_txi_deseq <- toRet3
      saveRDS(
        res_txi_deseq,
        file = file.path(new_dir,'res_txi_deseq.RDS')
      )
      openxlsx::write.xlsx(res_txi_deseq,
                           file = file.path(new_dir,'DEGs_full.xlsx'))

      # recalculating res_DEGs_txi_deseq
      toRet <- res_txi_deseq
      if(class(toRet)=="data.frame"){
        toRet <- toRet[which(toRet$padj<0.05),]
        if(nrow(toRet)==0) return({as.data.frame("oops, none significant")})
        toRet <- toRet[order(toRet$padj),]
        openxlsx::write.xlsx(as.data.frame(toRet), file = file.path(new_dir,'DEGs.xlsx'))
      } else {
        for (i in 1:length(toRet)){
          toRet[[i]] <- toRet[[i]][which(toRet[[i]]$padj<0.05),]
          if(nrow(toRet[[i]])==0){
            toRet[[i]] <- "oops, none significant"
          }
          else {
            toRet[[i]] <- toRet[[i]][order(toRet[[i]]$padj),]
          }
        }
        openxlsx::write.xlsx(toRet, file = file.path(new_dir,'DEGs.xlsx'))
      }
      res_DEGs_txi_deseq <- toRet
      saveRDS(
        res_DEGs_txi_deseq,
        file = file.path(new_dir,'res_DEGs_txi_deseq.RDS')
      )

      # recalculating GO_result
      toRet <- toRet2 <- toRet3 <- toWrite <- tytl <- list()
      if(class(res_DEGs_txi_deseq)=='list'){
        # multiple groups
        for (i in 1:length(res_DEGs_txi_deseq)){
          tytl[[i]] <- paste(strsplit(names(res_txi_deseq[i]), split = "_")[[1]][2:4], collapse=' ')
          if(is.null(dim(res_DEGs_txi_deseq[[i]]))){
            toRet[[i]] <- NULL
            toWrite[[i]] <- NULL
          } else {
            toRet[[i]] <- clusterProfiler::enrichGO(gene = res_DEGs_txi_deseq[[i]][,2],
                                                    keyType = "SYMBOL",
                                                    OrgDb = rv$OrgDeeBee,
                                                    ont = "BP",
                                                    pAdjustMethod = "BH",
                                                    qvalueCutoff = 0.05,
                                                    readable = TRUE)

            toWrite[[i]] <- as.data.frame(toRet[[i]])
            toWrite[[i]]$geneID <- gsub(pattern='/', replacement=' ', toWrite[[i]]$geneID)
          }
        }
        names(toRet) <- unlist(tytl)
        saveRDS(toRet, file = file.path(new_dir,'GO_result.RDS'))
        GO_result <- toRet
        names(toRet) <- gsub(pattern = 'Group ', replacement = '', x = names(toRet))
        names(toRet) <- substr(names(toRet), start = 1, stop = 30)
        openxlsx::write.xlsx(toWrite, file = file.path(new_dir,'GOs.xlsx'))
        return(toRet)

      } else {
        # 1 vs 1
        toRet <- clusterProfiler::enrichGO(gene = res_DEGs_txi_deseq[,2],
                                           keyType = "SYMBOL",
                                           OrgDb = rv$OrgDeeBee,
                                           ont = "BP",
                                           pAdjustMethod = "BH",
                                           qvalueCutoff = 0.05,
                                           readable = TRUE)
        saveRDS(toRet, file = file.path(new_dir,'GO_result.RDS'))
        GO_result <- toRet
        toWrite <- as.data.frame(toRet)
        toWrite$geneID <- gsub(pattern='/', replacement=' ', toWrite$geneID)
        openxlsx::write.xlsx(toWrite, file = file.path(new_dir,'GOs.xlsx'))
        return(toRet)
      }

      saveRDS(
        DESeq2::vst(txi_deseq),
        file = file.path(new_dir, "vst_data.RDS")
      )
      message("colData rebased")

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

      saveRDS(txi_subset, file.path(new_dir, "txi.RDS"))
      saveRDS(colData_subset, file.path(new_dir, "colData.RDS"))
      saveRDS(rv$referenceGenomeChoice,
              file.path(new_dir, "referenceGenomeChoice.RDS"))


      #res <- format_df_numbers(txi_subset$abundance, digits = 3)
      #res <- data.frame(ENSEMBL=rownames(res), res)
      #annots <- AnnotationDbi::select(rv$OrgDeeBee, keys=rownames(res),
      #                                columns="SYMBOL", keytype="ENSEMBL")
      #result <- merge(annots, res, by.x="ENSEMBL", by.y="ENSEMBL")
      #result <- rbind(c('-','-',colData_subset$Group),result)

      #openxlsx::write.xlsx(result, file = file.path(new_dir, 'TPMs.xlsx'))
      #saveRDS(
      #  result,
      #  file = file.path(new_dir, 'txi_tpms.RDS')
      #)

      reff <- as.character(colData_subset[which(colData_subset[,3]),2][1])
      if(length(which(colData_subset[,3]))>0) colData_subset$Group <- relevel(x = colData_subset$Group, ref = reff)
      txi_deseq <- DESeq2::DESeqDataSetFromTximport(txi_subset, colData = colData_subset, design = ~Group)

      saveRDS(
        txi_deseq,
        file = file.path(new_dir, "txi_deseq.RDS")
      )

      matr <- txi_deseq
      smallestGroupSize <- floor(min(table(colData_subset$Group)))
      keep <- rowSums(DESeq2::counts(matr) >= 10) >= smallestGroupSize
      matr <- matr[keep,]
      toRet <- DESeq2::DESeq(matr)
      txi_deseq_deseq <- toRet

      saveRDS(txi_deseq_deseq,
              file = file.path(new_dir, 'txi_deseq_deseq.RDS'))

      toRet <- toRet2 <- toRet3 <- toWrite <- tytl <- list()
      if(length(DESeq2::resultsNames(txi_deseq_deseq))<=2){
        toRet <- DESeq2::results(txi_deseq_deseq)
        toRet$log10padj <- log10(toRet$padj)
        toRet$Significance <- 'Non-significant'
        toRet$Significance[toRet$padj<0.05] <- 'Significant (padj < 0.05)'

        toRet2 <- data.frame(ENSEMBL=rownames(toRet),
                             toRet[,1:2],
                             abs_log2FoldChange=abs(toRet[,2]),toRet[,3:ncol(toRet)])
        annots <- AnnotationDbi::select(rv$OrgDeeBee, keys=rownames(toRet2),
                                        columns="SYMBOL", keytype="ENSEMBL")

        toRet3 <- merge(annots, toRet2, by.x="ENSEMBL", by.y="ENSEMBL")
        toRet3 <- toRet3[order(toRet3$padj),]

      }  else {
        toRet3 <- toRet <- list()
        for (i in 2:length(DESeq2::resultsNames(txi_deseq_deseq))){
          toRet[[(i-1)]] <- DESeq2::results(txi_deseq_deseq, name = DESeq2::resultsNames(txi_deseq_deseq)[i] )
          toRet[[(i-1)]]$log10padj <- log10( toRet[[(i-1)]]$padj)
          toRet[[(i-1)]]$Significance <- 'Non-significant'
          toRet[[(i-1)]]$Significance[ toRet[[(i-1)]]$padj<0.05] <- 'Significant (padj < 0.05)'

          toRet2 <- data.frame(ENSEMBL=rownames(toRet[[(i-1)]]),
                               toRet[[(i-1)]][,1:2],
                               abs_log2FoldChange=abs(toRet[[(i-1)]][,2]),toRet[[(i-1)]][,3:ncol(toRet[[(i-1)]])])

          annots <- AnnotationDbi::select(rv$OrgDeeBee, keys=rownames(toRet2),
                                          columns="SYMBOL", keytype="ENSEMBL")

          toRet3[[(i-1)]] <- merge(annots,
                                   toRet2,
                                   by.x="ENSEMBL",
                                   by.y="ENSEMBL")
          toRet3[[(i-1)]] <- toRet3[[(i-1)]][order(toRet3[[(i-1)]]$padj),]

        }
        names(toRet3) <- DESeq2::resultsNames(txi_deseq_deseq)[2:length(DESeq2::resultsNames(txi_deseq_deseq))]
      }
      res_txi_deseq <- toRet3
      saveRDS(
        res_txi_deseq,
        file = file.path(new_dir,'res_txi_deseq.RDS')
      )
      openxlsx::write.xlsx(res_txi_deseq,
                           file = file.path(new_dir,'DEGs_full.xlsx'))

      # recalculating res_DEGs_txi_deseq
      toRet <- res_txi_deseq
      if(class(toRet)=="data.frame"){
        toRet <- toRet[which(toRet$padj<0.05),]
        if(nrow(toRet)==0) return({as.data.frame("oops, none significant")})
        toRet <- toRet[order(toRet$padj),]
        openxlsx::write.xlsx(as.data.frame(toRet), file = file.path(new_dir,'DEGs.xlsx'))
      } else {
        for (i in 1:length(toRet)){
          toRet[[i]] <- toRet[[i]][which(toRet[[i]]$padj<0.05),]
          if(nrow(toRet[[i]])==0){
            toRet[[i]] <- "oops, none significant"
          }
          else {
            toRet[[i]] <- toRet[[i]][order(toRet[[i]]$padj),]
          }
        }
        openxlsx::write.xlsx(toRet, file = file.path(new_dir,'DEGs.xlsx'))
      }
      res_DEGs_txi_deseq <- toRet
      saveRDS(
        res_DEGs_txi_deseq,
        file = file.path(new_dir,'res_DEGs_txi_deseq.RDS')
      )

      # recalculating GO_result
      toRet <- toRet2 <- toRet3 <- toWrite <- tytl <- list()
      if(class(res_DEGs_txi_deseq)=='list'){
        # multiple groups
        for (i in 1:length(res_DEGs_txi_deseq)){
          tytl[[i]] <- paste(strsplit(names(res_txi_deseq[i]), split = "_")[[1]][2:4], collapse=' ')
          if(is.null(dim(res_DEGs_txi_deseq[[i]]))){
            toRet[[i]] <- NULL
            toWrite[[i]] <- NULL
          } else {
            toRet[[i]] <- clusterProfiler::enrichGO(gene = res_DEGs_txi_deseq[[i]][,2],
                                   keyType = "SYMBOL",
                                   OrgDb = rv$OrgDeeBee,
                                   ont = "BP",
                                   pAdjustMethod = "BH",
                                   qvalueCutoff = 0.05,
                                   readable = TRUE)

            toWrite[[i]] <- as.data.frame(toRet[[i]])
            toWrite[[i]]$geneID <- gsub(pattern='/', replacement=' ', toWrite[[i]]$geneID)
          }
        }
        names(toRet) <- unlist(tytl)
        saveRDS(toRet, file = file.path(new_dir,'GO_result.RDS'))
        GO_result <- toRet
        names(toRet) <- gsub(pattern = 'Group ', replacement = '', x = names(toRet))
        names(toRet) <- substr(names(toRet), start = 1, stop = 30)
        openxlsx::write.xlsx(toWrite, file = file.path(new_dir,'GOs.xlsx'))
        return(toRet)

      } else {
        # 1 vs 1
        toRet <- clusterProfiler::enrichGO(gene = res_DEGs_txi_deseq[,2],
                          keyType = "SYMBOL",
                          OrgDb = rv$OrgDeeBee,
                          ont = "BP",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05,
                          readable = TRUE)
        saveRDS(toRet, file = file.path(new_dir,'GO_result.RDS'))
        GO_result <- toRet
        toWrite <- as.data.frame(toRet)
        toWrite$geneID <- gsub(pattern='/', replacement=' ', toWrite$geneID)
        openxlsx::write.xlsx(toWrite, file = file.path(new_dir,'GOs.xlsx'))
        return(toRet)
      }

      saveRDS(
        DESeq2::vst(txi_deseq),
        file = file.path(new_dir, "vst_data.RDS")
      )
      message("colData rebased")



      output$forkingfeedback <- renderText(
        paste0("Successfully copied samples into folder ", foldr)
      )

    })

  })
}
