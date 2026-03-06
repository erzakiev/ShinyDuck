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
