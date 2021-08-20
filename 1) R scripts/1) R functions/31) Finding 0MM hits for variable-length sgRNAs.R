### 20th July 2021 ###




# Define functions --------------------------------------------------------

RowBindToDf <- function(df_list) {
  do.call(rbind.data.frame,
          c(df_list,
            list(stringsAsFactors = FALSE,
                 make.row.names = FALSE
            )))
}


FindAllSubRanges <- function(whole_length, min_length) {
  if (whole_length < min_length) {
    stop("'whole_length' must be at least as large as 'min_length'!")
  }
  results_list <- list(c(1, whole_length))
  if (min_length < whole_length) {
    smaller_lengths <- seq(whole_length - 1L, min_length, by = -1L)
    smaller_list <- lapply(smaller_lengths, function(x) {
      start_seq <- seq_len(whole_length - x + 1L)
      end_seq <- start_seq + x - 1L
      mapply(function(x, y) c(x, y), start_seq, end_seq, SIMPLIFY = FALSE)
    })
    results_list <- c(results_list, unlist(smaller_list, recursive = FALSE))
  }
  results_mat <- do.call(rbind, results_list)
  colnames(results_mat) <- c("Start", "End")
  results_mat <- cbind("Length" = results_mat[, "End"] - results_mat[, "Start"] + 1L,
                       results_mat
                       )
  return(results_mat)
}



GetAllSubSequences <- function(sequences_vec, min_length = 18L) {
  stopifnot(!(any(duplicated(sequences_vec))))
  sequence_lengths <- nchar(sequences_vec)
  unique_lengths <- sort(unique(sequence_lengths))
  if (any(unique_lengths < min_length)) {
    stop(paste0("Some of the sequences were shorter than 'min_length'",
                "(i.e. ", min_length, ")!"
                )
         )
  }
  unique_lengths_list <- lapply(unique_lengths, function(x) {
    sub_range_mat <- FindAllSubRanges(x, min_length)
    are_this_length <- sequence_lengths == x
    ranges_list <- lapply(seq_len(nrow(sub_range_mat)), function(y) {
      truncated_vec <- substr(sequences_vec[are_this_length],
                              sub_range_mat[y, "Start"],
                              sub_range_mat[y, "End"]
                              )
      results_df <- data.frame(
        "Original"  = sequences_vec[are_this_length],
        "Truncated" = truncated_vec,
        sub_range_mat[rep(y, sum(are_this_length)), , drop = FALSE],
        stringsAsFactors = FALSE
      )
    })
    ranges_df <- RowBindToDf(ranges_list)
  })
  unique_lengths_df <- RowBindToDf(unique_lengths_list)
  new_order <- order(match(unique_lengths_df[["Original"]], sequences_vec))
  unique_lengths_df <- unique_lengths_df[new_order, ]
  row.names(unique_lengths_df) <- NULL
  return(unique_lengths_df)
}




SearchFor0MMHits <- function(sequences_vec) {
  unique_sequences <- unique(sequences_vec)
  sequence_lengths <- nchar(unique_sequences)
  unique_lengths <- sort(unique(sequence_lengths))
  sequences_df_list <- lapply(unique_lengths, function(x) {
    message(paste0("Looking for matches for ", x, "-nucleotide sequences..."))
    results_df <- FindSequences(unique_sequences[sequence_lengths == x], max.mismatch = 0)
    message("\n")
    return(results_df)
  })
  results_df <- RowBindToDf(sequences_df_list)
  results_df <- results_df[, !(names(results_df) %in% c("Reference", "Num_MM"))]
  return(results_df)
}



FindPerfectHits <- function(input_sequences,
                            min_length = 19L,
                            valid_PAMs = "GG" #"GA", "AG"
                            ) {

  sub_sequences_df <- GetAllSubSequences(toupper(input_sequences), min_length = min_length)

  hits_df <- SearchFor0MMHits(sub_sequences_df[["Truncated"]])
  hits_df[["PAM"]] <- GetNGGPAM(hits_df)

  are_valid <- substr(hits_df[["PAM"]], 2, 3) %in% valid_PAMs

  hits_df <- hits_df[are_valid, ]
  row.names(hits_df) <- NULL

  expanded_hits_list <- lapply(seq_len(nrow(sub_sequences_df)), function(x) {
    this_sequence <- sub_sequences_df[["Truncated"]][[x]]
    are_this_sequence <- hits_df[["Sequence"]] == this_sequence
    if (any(are_this_sequence)) {
      results_df <- data.frame("Original"  = sub_sequences_df[["Original"]][[x]],
                               "Truncated" = this_sequence,
                               hits_df[are_this_sequence, names(hits_df) != "Sequence"],
                               stringsAsFactors = FALSE,
                               row.names = NULL
                               )
    } else {
      return(NULL)
    }
  })

  expanded_hits_df <- RowBindToDf(expanded_hits_list)
  expanded_hits_splits <- split(expanded_hits_df,
                                factor(expanded_hits_df[["Original"]],
                                       levels = unique(expanded_hits_df[["Original"]])
                                       )
                                )

  filtered_hits_splits <- lapply(expanded_hits_splits, function(x) {
    indices_vec <- seq_len(nrow(x))
    are_to_exclude <- vapply(indices_vec, function(y) {
      other_indices <- indices_vec[indices_vec != y]
      any(((x[["Start"]][other_indices] <  x[["Start"]][[y]]) & (x[["End"]][other_indices] >= x[["End"]][[y]])) |
          ((x[["Start"]][other_indices] <= x[["Start"]][[y]]) & (x[["End"]][other_indices] >  x[["End"]][[y]]))
          )
    }, logical(1))
    if (any(are_to_exclude)) {
      x <- x[!(are_to_exclude), ]
      row.names(x) <- NULL
    }
    return(x)
  })

  filtered_df <- RowBindToDf(filtered_hits_splits)
  filtered_df[["Location_string"]] <- MakeLocationStrings(filtered_df)
  sequence_fac <- factor(filtered_df[["Original"]],
                         levels = unique(filtered_df[["Original"]])
                         )

  hits_0MM_vec <- tapply(filtered_df[["Location_string"]], sequence_fac, function(x) {
    if (length(x) == 1) {
      return(x)
    } else {
      paste0(x, collapse = "; ")
    }
  })

  PAMs_vec <- tapply(filtered_df[["PAM"]], sequence_fac, function(x) {
    if (length(x) == 1) {
      return(x)
    } else {
      unique_PAMs <- unique(x)
      if (length(unique_PAMs) == 1) {
        return(unique_PAMs)
      } else {
        return(paste0(unique_PAMs, collapse = ", "))
      }
    }
  })

  summary_df <- data.frame(
    "Sequence" = levels(sequence_fac),
    "Num_0MM"  = tabulate(sequence_fac),
    "Loci_0MM" = hits_0MM_vec,
    "PAMs_0MM" = PAMs_vec,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  matches_vec <- match(toupper(input_sequences), summary_df[["Sequence"]])
  results_df <- summary_df[matches_vec, ]
  row.names(results_df) <- NULL
  results_df[["Sequence"]] <- input_sequences
  return(results_df)
}




Add0MMHits <- function(input_df, sequence_column, min_length = 19L) {

  unique_sequences <- unique(input_df[, sequence_column])
  loci_0MM_df <- FindPerfectHits(unique_sequences, min_length = min_length)

  matches_vec <- match(input_df[[sequence_column]],
                       loci_0MM_df[["Sequence"]]
                       )

  matched_df <- loci_0MM_df[matches_vec, names(loci_0MM_df) != "Sequence"]
  row.names(matched_df) <- NULL
  results_df <- data.frame(input_df,
                           matched_df,
                           stringsAsFactors = FALSE
                           )
  results_df[["Num_0MM"]][is.na(results_df[["Num_0MM"]])] <- 0L
  return(results_df)
}



