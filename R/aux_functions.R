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
  rownames(df) <- as.vector(convertVectEns2Symb(rownames(df), orgdb))
  return(df)
}


calculate_txi_tpm <- function(txi, orgdb, colData){

  res<-txi$abundance
  res <- data.frame(ENSEMBL=rownames(res), res)
  annots <- convertVectEns2Symb(rownames(res), orgdb = orgdb)
  result <- cbind(annots, res)
  result <- rbind(c('-','-',colData$data$Group),result)
  return(result)

}

res_counts<-txi()$counts
res_counts <- data.frame(ENSEMBL=rownames(res_counts), res_counts)
annots2 <- AnnotationDbi::select(OrgDeeBee(), keys=rownames(res_counts),
                                 columns="SYMBOL", keytype="ENSEMBL")
result_counts <- merge(annots2, res_counts, by.x="ENSEMBL", by.y="ENSEMBL")
result_counts <- rbind(c('-','-',colDatt$data$Group),result_counts)
openxlsx::write.xlsx(result_counts, file = paste0(ProjFolderFull(), '/Counts.xlsx'))


