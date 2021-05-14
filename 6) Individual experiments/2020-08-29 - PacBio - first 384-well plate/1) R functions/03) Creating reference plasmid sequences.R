### 20th October 2020 ###





# Import packages and source code -----------------------------------------

library("Biostrings")




# Define functions --------------------------------------------------------

ReplaceNNNLines <- function(lines_vec, sequences_vec) {

  ### Extract the positions of the 4 guides

  all_matches <- gregexpr("n+[n ]*", lines_vec)

  match_positions <- lapply(all_matches, function(x) as.integer(x))
  match_lengths <- lapply(all_matches, function(x) attributes(x)[["match.length"]])

  have_matches <- !(vapply(match_positions, function(x) identical(x, -1L), logical(1)))

  ### Group locations that span two lines

  match_numbers <- rep(NA_integer_, sum(have_matches))
  current_number <- 0L
  for (i in seq_along(match_numbers)) {
    if (i == 1) {
      previous_reaches_line_end <- FALSE
    } else {
      previous_reaches_line_end <- grepl("n+[n ]*$", lines_vec[have_matches][[i - 1]])
    }
    if (previous_reaches_line_end) {
      starts_from_line_start <- grepl("^[n 0-9]*n", lines_vec[have_matches][[i]])
      if (!(starts_from_line_start)) {
        current_number <- current_number + 1L
      }
    } else {
      current_number <- current_number + 1L
    }
    match_numbers[[i]] <- current_number
  }

  ### Prepare for replacements

  stopifnot(all(match_numbers %in% 1:6))

  indices_list <- split(seq_len(sum(have_matches)), match_numbers)
  new_lines_vec <- rep(NA_character_, sum(have_matches))

  matches_per_line <- vapply(indices_list,
                             function(x) max(lengths(match_positions[have_matches][x])),
                             integer(1)
                             )
  sequences_list <- split(sequences_vec, rep(seq_along(indices_list), matches_per_line))

  ### Replace the 4 guide sequences

  for (i in seq_along(indices_list)) {
    sum_replaced <- 0L
    for (j in indices_list[[i]]) {
      match_index <- which(have_matches)[[j]]
      current_line <- lines_vec[[match_index]]
      num_per_line <- length(match_positions[[match_index]])
      for (k in seq_len(num_per_line)) {
        nnn_start <- match_positions[[match_index]][[k]]
        nnn_length <- match_lengths[[match_index]][[k]]
        nnn_end <- nnn_start + nnn_length - 1L
        nnn_string <- substr(current_line, nnn_start, nnn_end)
        nnn_chars <- strsplit(nnn_string, "", fixed = TRUE)[[1]]
        replacement_chars <- rep(NA_character_, nnn_length)
        for (char_index in seq_len(nnn_length)) {
          if (nnn_chars[[char_index]] == " ") {
            replacement_chars[[char_index]] <- " "
          } else {
            sum_replaced <- sum_replaced + 1L
            replacement_chars[[char_index]] <- substr(sequences_list[[i]][[k]],
                                                      sum_replaced,
                                                      sum_replaced
                                                      )
          }
        }
        replacement_string <- paste0(replacement_chars, collapse = "")
        if (nnn_start != 1) {
          string_before <- substr(current_line, 1, nnn_start - 1L)
        } else {
          string_before <- ""
        }
        line_length <- nchar(current_line)
        if (nnn_end != line_length) {
          string_after <- substr(current_line, nnn_end + 1L, line_length)
        } else {
          string_after <- ""
        }
        current_line <- paste0(string_before, replacement_string, string_after)
        if ((num_per_line == 2) && (sum_replaced == nchar(sequences_list[[i]][[k]]))) {
          sum_replaced <- 0
        }
      }
      new_lines_vec[[j]] <- current_line
    }
  }

  results_vec <- lines_vec
  results_vec[have_matches] <- new_lines_vec
  return(results_vec)
}



ReadInPlasmid_gbk <- function(file_path) {
  lines_vec <- scan(file = file_path, what = character(), sep = "\n")

  first_index <- which(lines_vec == "ORIGIN") + 1L
  preceding_lines <- lines_vec[seq_len(first_index - 1)]
  last_index <- first_index + 37L
  following_seq <- seq(from = last_index + 1, to = length(lines_vec))
  following_lines <- lines_vec[following_seq]
  sequence_seq <- seq(from = first_index, to = last_index)
  sequence_lines <- lines_vec[sequence_seq]

  result_list <- list(
    "preceding" = preceding_lines,
    "sequence"  = sequence_lines,
    "following" = following_lines
  )
  return(result_list)
}






AddReferenceSequences <- function(sg_df, tracRNAs, promoters, plasmid) {

  stopifnot(length(tracRNAs) == 4)
  stopifnot(length(promoters) == 4)
  sg_columns <- paste0("Sequence_sg", 1:4)
  stopifnot(all(sg_columns %in% names(sg_df)))

  rev_tracRNAs <- as.character(reverseComplement(DNAStringSet(tracRNAs)))

  sg_tracRNA_mat <- do.call(cbind, lapply(1:4, function(x) {
    paste0(sg_df[[sg_columns[[x]]]], rev_tracRNAs[[x]], "TTTTT")
  }))
  colnames(sg_tracRNA_mat) <- paste0("sg_cr_", 1:4)

  sg_with_promoters_mat <- do.call(cbind, lapply(1:4, function(x) {
    paste0(toupper(promoters[[x]]), sg_tracRNA_mat[, x])
  }))
  colnames(sg_with_promoters_mat) <- paste0("sg_cr_pr_", 1:4)


  plasmids_vec <- vapply(seq_len(nrow(sg_df)), function(x) {
    sg_N <- paste0(rep("N", 20), collapse = "")
    sg_sequences <- vapply(1:4, function(y) {
      sg_df[[paste0("Sequence_sg", y)]][[x]]
    }, "")
    result_string <- sub(sg_N, sg_sequences[[1]], plasmid, fixed = TRUE)
    result_string <- sub(sg_N, sg_sequences[[2]], result_string, fixed = TRUE)
    result_string <- sub(sg_N, sg_sequences[[3]], result_string, fixed = TRUE)
    result_string <- sub(sg_N, sg_sequences[[4]], result_string, fixed = TRUE)
    return(result_string)
  }, "")

  sg_sequences_df <- data.frame(
    sg_df,
    sg_tracRNA_mat,
    sg_with_promoters_mat,
    "Whole_plasmid" = plasmids_vec,
    stringsAsFactors = FALSE
  )
  return(sg_sequences_df)
}








