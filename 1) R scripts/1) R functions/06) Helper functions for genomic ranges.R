### 10 August 2019 ###





# Define functions --------------------------------------------------------

CheckRangesDf <- function(ranges_df) {
  required_columns <- c("Chromosome", "Strand", "Start", "End")
  if (!(all(required_columns %in% names(ranges_df)))) {
    stop(paste0("The data frame did not contain all required columns (", paste0(required_columns, collapse = ", "), ")"))
  }
}





