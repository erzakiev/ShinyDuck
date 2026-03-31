format_df_numbers <- function(df,
                              digits = 2,
                              sci_threshold = 1e-4) {

  num_cols <- sapply(df, function(x) is.numeric(x) && !is.integer(x))

  df[num_cols] <- lapply(
    df[num_cols],
    format_numbers,
    digits = digits,
    sci_threshold = sci_threshold
  )

  df[num_cols] <- lapply(
    df[num_cols],
    as.numeric
  )

  signif_cols <- grep(colnames(df), pattern = 'signif', value = T, ignore.case = T)
  if(length(signif_cols)>0){
    for (i in 1:length(signif_cols)){
      df[,signif_cols[i]] <- as.factor(df[,signif_cols[i]])
    }
  }
  df
}

format_numbers <- function(x,
                           digits = 2,
                           sci_threshold = 1e-4) {
  ifelse(
    abs(x) < sci_threshold & x != 0,
    formatC(x, format = "e", digits = digits),
    formatC(x, format = "f", digits = digits)
  )
}

colorize_df <- function(df){

}

strsplits <- function(x, splits = c(" ", ",", "/", "\n", "\t"))
{
  for (split in splits) {
    x <- unlist(strsplit(x, split, fixed = TRUE))
  }

  x <- trimws(x)
  x <- x[x != ""]

  unique(x)
}

convertVectEns2Symb <- function(vect, orgdb){
  annots <- AnnotationDbi::mapIds(
      orgdb,
      keys = vect,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
  )
  annots[is.na(annots)] <- vect[is.na(annots)]
  return(annots)
}

convertDfRownamesEns2Symb <- function(df, orgdb){
  rownames(df) <- make.unique(as.vector(convertVectEns2Symb(rownames(df), orgdb)))
  return(df)
}


#calculate_txi_tpm <- function(txi, colData, orgdb){
#
#  res<-txi$abundance
#  res <- data.frame(ENSEMBL=rownames(res), res)
#  annots <- convertVectEns2Symb(rownames(res), orgdb = orgdb)
#  result <- cbind(annots, res)
#  result <- rbind(c('-','-',colData$Group),result)
#  return(result)
#
#}

calculate_txi_tpm <- function(txi, orgdb){

  res<-txi$abundance
  res <- convertDfRownamesEns2Symb(res, orgdb)
  return(res)

}

calculate_res_txi_deseq <- function(txi_deseq_deseq, orgdb){

  toRet <- toRet2 <- toRet3 <- toWrite <- tytl <- list()
  if(length(DESeq2::resultsNames(txi_deseq_deseq))<=2){
    toRet <- DESeq2::results(txi_deseq_deseq)
    toRet$log10padj <- log10(toRet$padj)
    toRet$Significance <- 'Non-significant'
    toRet$Significance[toRet$padj<0.05] <- 'Significant (padj < 0.05)'

    toRet2 <- data.frame(ENSEMBL=rownames(toRet),
                         toRet[,1:2],
                         abs_log2FoldChange=abs(toRet[,2]),toRet[,3:ncol(toRet)])
    annots <- AnnotationDbi::select(orgdb, keys=rownames(toRet2),
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

      annots <- AnnotationDbi::select(orgdb, keys=rownames(toRet2),
                                      columns="SYMBOL", keytype="ENSEMBL")

      toRet3[[(i-1)]] <- merge(annots,
                               toRet2,
                               by.x="ENSEMBL",
                               by.y="ENSEMBL")
      toRet3[[(i-1)]] <- toRet3[[(i-1)]][order(toRet3[[(i-1)]]$padj),]

    }
    names(toRet3) <- DESeq2::resultsNames(txi_deseq_deseq)[2:length(DESeq2::resultsNames(txi_deseq_deseq))]
  }
  return(toRet3)
}

calculate_res_DEGs_txi_deseq <- function(res_txi_deseq){

  toRet <- res_txi_deseq
  if(class(toRet)=="data.frame"){
    toRet <- toRet[which(toRet$padj<0.05),]
    if(nrow(toRet)==0) return({as.data.frame("oops, none significant")})
    toRet <- toRet[order(toRet$padj),]
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
  }
  return(toRet)
}

calculate_GO_result <- function(res_DEGs_txi_deseq, orgdb){
  toRet <- toRet2 <- toRet3 <- toWrite <- tytl <- list()
  if(class(res_DEGs_txi_deseq)=='list'){
    # multiple groups
    for (i in 1:length(res_DEGs_txi_deseq)){
      tytl[[i]] <- gsub(pattern = 'Group_', replacement = '', names(res_DEGs_txi_deseq[i]))
      if(is.null(dim(res_DEGs_txi_deseq[[i]]))){
        toRet[[i]] <- NULL
        toWrite[[i]] <- NULL
      } else {
        incProgress(0.05, detail = 'Recalculating the GO enrichments')
        toRet[[i]] <- clusterProfiler::enrichGO(gene = res_DEGs_txi_deseq[[i]][,2],
                                                keyType = "SYMBOL",
                                                OrgDb = orgdb,
                                                ont = "BP",
                                                pAdjustMethod = "BH",
                                                qvalueCutoff = 0.05,
                                                readable = TRUE)
      }
    }
    names(toRet) <- unlist(tytl)
  } else {
    # 1 vs 1
    toRet <- clusterProfiler::enrichGO(gene = res_DEGs_txi_deseq[,2],
                                       keyType = "SYMBOL",
                                       OrgDb = orgdb,
                                       ont = "BP",
                                       pAdjustMethod = "BH",
                                       qvalueCutoff = 0.05,
                                       readable = TRUE)
  }
  return(toRet)
}

reshape_GO_result_for_xlsx <- function(res_DEGs_txi_deseq, orgdb){
  toRet <- toRet2 <- toRet3 <- toWrite <- tytl <- list()
  if(class(res_DEGs_txi_deseq)=='list'){
    # multiple groups
    for (i in 1:length(res_DEGs_txi_deseq)){
      tytl[[i]] <- gsub(pattern = 'Group_', replacement = '', names(res_DEGs_txi_deseq[i]))
      if(is.null(dim(res_DEGs_txi_deseq[[i]]))){
        toRet[[i]] <- NULL
        toWrite[[i]] <- NULL
      } else {
        incProgress(0.05, detail = 'Recalculating the GO enrichments')
        toRet[[i]] <- clusterProfiler::enrichGO(gene = res_DEGs_txi_deseq[[i]][,2],
                                                keyType = "SYMBOL",
                                                OrgDb = orgdb,
                                                ont = "BP",
                                                pAdjustMethod = "BH",
                                                qvalueCutoff = 0.05,
                                                readable = TRUE)

        toWrite[[i]] <- as.data.frame(toRet[[i]])
        toWrite[[i]]$geneID <- gsub(pattern='/', replacement=' ', toWrite[[i]]$geneID)
      }
    }
    names(toRet) <- unlist(tytl)
    return(toRet)

    incProgress(0.05, detail = 'Saving the GO enrichments')
    saveRDS(toRet, file = file.path(new_dir,'GO_result.RDS'))
    names(toRet) <- gsub(pattern = 'Group ', replacement = '', x = names(toRet))
    names(toRet) <- substr(names(toRet), start = 1, stop = 30)
    incProgress(0.05, detail = 'Writing the GOs.xlsx')
    openxlsx::write.xlsx(toWrite, file = file.path(new_dir,'GOs.xlsx'))
    return(toRet)

  } else {
    # 1 vs 1
    incProgress(0.05, detail = 'Recalculating the GO enrichments')
    toRet <- clusterProfiler::enrichGO(gene = res_DEGs_txi_deseq[,2],
                                       keyType = "SYMBOL",
                                       OrgDb = orgdb,
                                       ont = "BP",
                                       pAdjustMethod = "BH",
                                       qvalueCutoff = 0.05,
                                       readable = TRUE)
    incProgress(0.05, detail = 'Saving the GO enrichments')
    saveRDS(toRet, file = file.path(new_dir,'GO_result.RDS'))
    toWrite <- as.data.frame(toRet)
    toWrite$geneID <- gsub(pattern='/', replacement=' ', toWrite$geneID)
    openxlsx::write.xlsx(toWrite, file = file.path(new_dir,'GOs.xlsx'))
  }
  return(toWrite)
}

calculate_txi_deseq_deseq <- function(txi_deseq, colData){
  matr <- txi_deseq
  smallestGroupSize <- floor(min(table(colData$Group)))
  keep <- rowSums(DESeq2::counts(matr) >= 10) >= smallestGroupSize
  matr <- matr[keep,]
  toRet <- DESeq2::DESeq(matr)
  return(toRet)
}


color_numeric_column <- function(dt, data, column,
                                 palette=c("#4169E1","white","red")){

  rng <- range(data[[column]], na.rm=TRUE)

  dt |>
    DT::formatStyle(
      column,

      backgroundColor = DT::styleInterval(
        seq(rng[1], rng[2], length.out = 999),

        colorRampPalette(palette)(1000)
      )
    )
}

color_tpm_table <- function(
    dt,
    data,
    global_scale = TRUE,
    log_scale = TRUE,
    palette = c("white","#fff7bc","#fec44f","#d95f0e","#993404"),
    bins = 20
){

  numeric_cols <- names(data)[sapply(data, is.numeric)]

  mat <- data[numeric_cols]

  if(log_scale){
    mat <- log10(mat + 1)
  }

  if(global_scale){

    all_values <- unlist(mat)

    rng <- range(all_values, na.rm = TRUE)

    breaks <- seq(rng[1], rng[2], length.out = bins)

    colors <- colorRampPalette(palette)(bins + 1)

    for(col in numeric_cols){

      dt <- dt |>
        DT::formatStyle(
          col,

          backgroundColor = DT::styleInterval(
            breaks,
            colors
          )
        )
    }

  } else {

    for(col in numeric_cols){

      rng <- range(mat[[col]], na.rm = TRUE)

      breaks <- seq(rng[1], rng[2], length.out = bins)

      colors <- colorRampPalette(palette)(bins + 1)

      dt <- dt |>
        DT::formatStyle(
          col,

          backgroundColor = DT::styleInterval(
            breaks,
            colors
          )
        )
    }

  }

  dt
}

make_TPM_container <- function(df, colData){

  #saveRDS(df, '~/df.RDS')
  sample_names <- colnames(df)[-1]
  #saveRDS(sample_names, '~/sample_names.RDS')
  #saveRDS(colData, '~/colData.RDS')

  sample_groups <- as.character(colData$Group[
    match(sample_names, colData$Sample)
  ])
  #saveRDS(sample_groups, '~/sample_groups.RDS')

  group_lengths <- rle(sample_groups)

  htmltools::withTags(

    table(

      thead(

        tr(

          th(colspan = 1,
             ""), # name

          lapply(

            seq_along(group_lengths$values),

            function(i){

              th(
                colspan = group_lengths$lengths[i],

                style = paste0(
                  "text-align:center;",
                  if(i < length(group_lengths$values))
                    "border-right:2px solid black;"
                  else
                    ""
                ),

                group_lengths$values[i]
              )

            }

          )

        ),

        tr(
          lapply(colnames(df), th)

        )

      )

    )

  )

}


add_group_separators <- function(dt, df, colData){

  sample_names <- colnames(df)[2:ncol(df)]

  sample_groups <- colData[
    match(sample_names, colData$Sample),
    "Group",
    drop = TRUE
  ]

  sample_groups <- as.character(sample_groups)

  r <- rle(sample_groups)

  group_end_positions <- cumsum(r$lengths)

  group_end_positions <- group_end_positions[-length(group_end_positions)]

  for(pos in group_end_positions){

    colname <- colnames(df)[-1][pos]

    dt <- dt |>
      DT::formatStyle(
        columns = colname,
        borderRight = "2px solid black"
      )

  }

  dt
}
