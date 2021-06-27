### 15th June 2021 ###



# Define functions --------------------------------------------------------

AddNumOccurrences <- function(sg_df) {

  all_sequences <- unlist(sg_df[, paste0("Sequence_sg", 1:4)])
  num_occurrences_mat <- do.call(cbind,
                                 lapply(sg_df[, paste0("Sequence_sg", 1:4)],
                                        function(x) vapply(x, function(y) sum(y == all_sequences), integer(1), USE.NAMES = FALSE)
                                 ))
  colnames(num_occurrences_mat) <- paste0("Num_occurrences_sg", 1:4)

  results_df <- data.frame(sg_df, num_occurrences_mat, stringsAsFactors = FALSE)
  return(results_df)
}
