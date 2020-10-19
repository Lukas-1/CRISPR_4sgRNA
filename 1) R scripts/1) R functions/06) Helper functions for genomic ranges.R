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



LocationsDfToGPosObject <- function(locations_df) {
  GPos_object <- GPos(
    seqnames = locations_df[["Chromosome"]],
    pos      = locations_df[["Cut_location"]],
    strand   = locations_df[["Strand"]]
  )
  return(GPos_object)
}





LocationStringToDf <- function(location_char_vec) {

  first_split_list <- strsplit(location_char_vec, ":", fixed = TRUE)
  first_half_vec <- sapply(first_split_list, "[[", 1)
  second_half_vec <- sapply(first_split_list, "[[", 2)

  nchar_first_half_vec <- nchar(first_half_vec)
  strand_pos <- nchar_first_half_vec - 1
  location_splits <- strsplit(second_half_vec, "-", fixed = TRUE)

  results_df <- data.frame(
    "Chromosome" = substr(first_half_vec, 1, nchar_first_half_vec - 3),
    "Strand"     = substr(first_half_vec, strand_pos, strand_pos),
    "Start"      = as.integer(sapply(location_splits, "[[", 1)),
    "End"        = as.integer(sapply(location_splits, "[[", 2)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  return(results_df)
}



TruncateLongEntries <- function(char_vec, use_sep = "; ", max_length = 10L) {
  splits_list <- strsplit(char_vec, use_sep, fixed = TRUE)
  TruncateLongEntriesSplits(splits_list, use_sep = use_sep, max_length = max_length, char_vec = char_vec)
}


TruncateLongEntriesSplits <- function(splits_list, use_sep = "; ", max_length = 10L, char_vec = NULL) {
  are_too_long <- lengths(splits_list) > max_length
  if (is.null(char_vec)) {
    char_vec <- rep(NA_character_, length(splits_list))
    char_vec[!(are_too_long)] <- vapply(splits_list[!(are_too_long)],
                                        function(x) paste0(x, collapse = use_sep),
                                        ""
                                        )
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





