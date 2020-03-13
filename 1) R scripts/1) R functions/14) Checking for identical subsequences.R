### 2nd December 2019 ###




# Unrelated helper functions ----------------------------------------------

CheckThatFactorIsInOrder <- function(my_factor) {
  stopifnot(identical(length(unique(my_factor)),
                      length(rle(as.integer(my_factor))[["lengths"]])
                      )
            )
}





# Functions for finding identical subsequences ----------------------------

SplitIntoSubstrings <- function(string, substring_length = 8L) {
  if (nchar(string) <= substring_length) {
    return(string)
  } else {
    start_positions <- seq_len(nchar(string) - substring_length)
    stop_positions <- seq(substring_length, nchar(string))
    sub_strings <- substr(rep.int(string, length(start_positions)), start_positions, stop_positions)
    return(sub_strings)
  }
}


NumHomologousPairs <- function(char_vec, substring_length = 8L) {
  comb_mat <- combn(length(char_vec), 2)
  comb_list <- lapply(seq_len(ncol(comb_mat)), function(x) comb_mat[, x])
  substrings_list <- lapply(char_vec, SplitIntoSubstrings, substring_length = substring_length)
  have_homologies <- vapply(comb_list, function(x) any(substrings_list[[x[[1]]]] %in% substrings_list[[x[[2]]]]), logical(1))
  return(sum(have_homologies))
}


ReturnIdenticalSubsequences <- function(char_vec, substring_length = 8L) {
  comb_mat <- combn(length(char_vec), 2)
  comb_list <- lapply(seq_len(ncol(comb_mat)), function(x) comb_mat[, x])
  substrings_list <- lapply(char_vec, SplitIntoSubstrings, substring_length = substring_length)
  homologies_list <- lapply(comb_list, function(x) intersect(substrings_list[[x[[1]]]], substrings_list[[x[[2]]]]))
  results_vec <- unlist(homologies_list)
  if (length(results_vec) == 0) {
    results_vec <- NA
  }
  return(results_vec)
}



LongestSharedSubsequence <- function(char_vec) {
  if (length(char_vec) < 2) {
    return(NA_integer_)
  }
  char_vec <- toupper(char_vec)
  max_num <- min(vapply(char_vec, nchar, integer(1)))
  nums_shared <- seq(from = max_num, to = 1L, by = -1L)
  for (num_shared in nums_shared) {
    num_homologies <- NumHomologousPairs(char_vec, num_shared)
    if (num_homologies > 0) {
      break
    }
  }
  return(num_shared)
}

VectorizedLongestSubsequenceSingle <- function(char_vec) {
  use_seq <- seq_along(char_vec)
  results_vec <- vapply(use_seq, function(x) LongestSubsequenceSingle(char_vec[[x]], char_vec[use_seq[use_seq != x]]), integer(1))
  return(results_vec)
}

LongestSubsequenceSingle <- function(single_char, char_vec) {
  assign("delete_single_char", single_char, envir = globalenv())
  stopifnot(length(single_char) == 1)
  assign("delete_char_vec", char_vec, envir = globalenv())
  longest_vec <- vapply(char_vec, function(x) LongestSharedSubsequence(c(single_char, x)), integer(1))
  return(max(longest_vec))
}



CheckForIdenticalSubsequences <- function(CRISPR_df, substring_length = 8L, ID_column = "AltTSS_ID") {
  are_not_controls <- CRISPR_df[["Is_control"]] %in% "No"
  are_top_four <- CRISPR_df[["Rank"]] %in% 1:4
  unique_IDs <- unique(CRISPR_df[[ID_column]][are_not_controls & are_top_four])
  are_homologous <- rep.int(FALSE, nrow(CRISPR_df))
  for (ID in unique_IDs) {
    are_this_ID <- (CRISPR_df[[ID_column]] == ID) & are_top_four
    if (sum(are_this_ID) > 1) {
      is_homologous <- NumHomologousPairs(CRISPR_df[["sgRNA_sequence"]][are_this_ID], substring_length = substring_length) > 0
      if (is_homologous) {
        are_homologous[are_this_ID] <- TRUE
      }
    }
  }
  return(are_homologous)
}












