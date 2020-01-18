### 15th January 2020 ###



# Define functions --------------------------------------------------------

AppendIDsWithoutEntrezs <- function(entrez_IDs_list, CRISPR_df) {
  have_no_entrez <- is.na(CRISPR_df[, "Entrez_ID"]) &
                    (CRISPR_df[, "Is_control"] == "No")
  IDs_list <- entrez_IDs_list
  IDs_list[[length(IDs_list)]] <- c(IDs_list[[length(IDs_list)]],
                                    unique(CRISPR_df[have_no_entrez, "Combined_ID"])
                                    )
  return(IDs_list)
}


CombineDfChunks <- function(df_list, max_num_per_chunk = 12000L) {
  current_sum <- nrow(df_list[[1]])
  current_chunk <- 1L
  chunk_vec <- c(1L, rep(NA_integer_, length(df_list) - 1))
  for (i in seq(from = 2, to = length(df_list), by = 1)) {
    current_nrow <- nrow(df_list[[i]])
    if ((current_sum + current_nrow) > max_num_per_chunk) {
      current_chunk <- current_chunk + 1L
      current_sum <- current_nrow
    } else {
      current_sum <- current_sum + current_nrow
    }
    chunk_vec[[i]] <- current_chunk
  }
  result_df_list <- tapply(df_list, chunk_vec, function(x) do.call(rbind.data.frame, c(x, list(stringsAsFactors = FALSE, make.row.names = FALSE))))
  names(result_df_list) <- tapply(sub("_TF", "", names(df_list)), chunk_vec, paste0, collapse = "")
  return(result_df_list)
}