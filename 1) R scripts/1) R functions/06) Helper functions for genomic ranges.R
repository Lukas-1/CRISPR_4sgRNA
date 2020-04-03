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
    ranges_df[["Chromosome"]] <- ifelse(ranges_df[["Chromosome"]] == "M", "MT", ranges_df[["Chromosome"]])
  }
  GRanges_object <- GRanges(
    seqnames = ranges_df[["Chromosome"]],
    ranges   = IRanges(start = ranges_df[["Start"]], end = ranges_df[["End"]]),
    strand   = ranges_df[["Strand"]]
  )
  return(GRanges_object)
}


TruncateLongEntries <- function(char_vec, use_sep = "; ", max_length = 10L) {
  splits_list <- strsplit(char_vec, use_sep, fixed = TRUE)
  TruncateLongEntriesSplits(splits_list, use_sep = use_sep, max_length = max_length, char_vec = char_vec)
}


TruncateLongEntriesSplits <- function(splits_list, use_sep = "; ", max_length = 10L, char_vec = NULL) {
  are_too_long <- lengths(splits_list) > max_length
  if (is.null(char_vec)) {
    char_vec <- rep(NA_character_, length(splits_list))
    char_vec[!(are_too_long)] <- vapply(splits_list[!(are_too_long)], function(x) paste0(x, collapse = use_sep), "")
  }
  char_vec[are_too_long] <- vapply(splits_list[are_too_long], function(x) {
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





