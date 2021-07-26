### 15th June 2021 ###



# Define functions --------------------------------------------------------

AddNumOccurrences <- function(sg_df) {
  all_sequences <- unlist(sg_df[, paste0("Sequence_sg", 1:4)])
  all_sequences_table <- table(all_sequences)
  num_occurrences_mat <- apply(sg_df[, paste0("Sequence_sg", 1:4)], 2, function(x) {
    matches_vec <- match(x, names(all_sequences_table))
    as.integer(all_sequences_table[matches_vec])
  })
  colnames(num_occurrences_mat) <- paste0("Num_occurrences_sg", 1:4)
  results_df <- data.frame(sg_df, num_occurrences_mat, stringsAsFactors = FALSE)
  return(results_df)
}
