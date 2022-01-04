### 29th December 2021 ###



# Define functions --------------------------------------------------------

AddNewZMWs <- function(use_ccs_df, use_align_df, pool_number) {

  are_this_pool <- use_ccs_df[, "Pool"] %in% pool_number
  ccs_read_IDs <- paste0(use_ccs_df[["Combined_ID"]][are_this_pool], "__",
                         use_ccs_df[["Original_ZMW"]][are_this_pool]
                         )

  align_read_IDs <- paste0(use_align_df[["Combined_ID"]], "__",
                           use_align_df[["Original_ZMW"]]
                           )

  are_included <- align_read_IDs %in% ccs_read_IDs

  matches_vec <- match(align_read_IDs[are_included], ccs_read_IDs)
  stopifnot(!(anyNA(matches_vec)))

  ccs_columns <- c("Run", "Pool", "ZMW", "Original_ZMW")

  results_df <- data.frame(
    use_ccs_df[are_this_pool, ][matches_vec, ccs_columns],
    use_align_df[are_included, !(names(use_align_df) %in% ccs_columns)],
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  return(results_df)
}
