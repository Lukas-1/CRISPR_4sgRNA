### 10 August 2019 ###





# Define functions --------------------------------------------------------

CheckRangesDf <- function(ranges_df) {
  required_columns <- c("Chromosome", "Strand", "Start", "End")
  if (!(all(required_columns %in% names(ranges_df)))) {
    stop(paste0("The data frame did not contain all required columns (", paste0(required_columns, collapse = ", "), ")"))
  }
}


RangesDfToGRangesObject <- function(ranges_df, strip_Chr = FALSE) {
  if (strip_Chr) {
    ranges_df[["Chromosome"]] <- sub("chr", "", ranges_df[["Chromosome"]], fixed = TRUE)
  }
  GRanges_object <- GRanges(
    seqnames = ranges_df[["Chromosome"]],
    ranges   = IRanges(start = ranges_df[["Start"]], end = ranges_df[["End"]]),
    strand   = ranges_df[["Strand"]]
  )
  return(GRanges_object)
}


TruncateLongEntries <- function(char_vec, use_sep = "; ", max_length = 10L) {
  splits <- strsplit(char_vec, use_sep, fixed = TRUE)
  are_too_long <- lengths(splits) > max_length
  char_vec[are_too_long] <- vapply(splits[are_too_long], function(x) {
    all_the_same <- length(unique(x)) == 1
    result <- paste0("Truncated (", length(x), " entries)...")
    if (all_the_same) {
      result <- paste0(result, " All entries are: ", x[[1]])
    } else {
      first_3_the_same <- length(unique(x[1:3])) == 1
      if (first_3_the_same) {
        result <- paste0(result, " The first 3 are all: ", x[[1]])
      } else {
        result <- paste0(result, " The first 3 are: ", paste0(x[1:3], collapse = use_sep))
      }
    }
    return(result)
  }, "")
  return(char_vec)
}
