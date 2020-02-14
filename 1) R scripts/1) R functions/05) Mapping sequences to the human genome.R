### 26 July 2019 ###





# Import packages and source code -----------------------------------------

library("BSgenome.Hsapiens.UCSC.hg38")




# Loading global variables ------------------------------------------------

chromosome_names <- paste0("chr", c(as.character(1:22), "X", "Y"))
message("Loading the human genome into RAM...")
chromosome_sequences_list <- sapply(chromosome_names, function(x) BSgenome.Hsapiens.UCSC.hg38[[x]], simplify = FALSE)





# Define functions --------------------------------------------------------

XStringViewsToDf <- function(matches_object, chromosome, strand, sequence) {
  if (strand == "-") {
    sequences_vec <- as.character(reverseComplement(DNAStringSet(as.character(matches_object))))
  } else {
    sequences_vec <- as.character(matches_object)
  }
  results_df <- data.frame(
    "Reference"  = sequence,
    "Sequence"   = sequences_vec,
    "Num_MM"     = ifelse(sequences_vec == sequence, 0L, 1L),
    "Chromosome" = chromosome,
    "Strand"     = strand,
    "Start"      = start(matches_object),
    "End"        = end(matches_object),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  return(results_df)
}



FindSequence <- function(my_sequence) {

  message("Looking for matches for ", my_sequence, " in the Hsapiens.UCSC.hg38 genome...")

  DS_sequence <- DNAString(my_sequence)
  DS_sequence_reversed <- reverseComplement(DS_sequence)
  results_list <- list()

  for (i in seq_along(chromosome_names)) {

    message(paste0("Searching chromosome ", substr(chromosome_names[[i]], 4, nchar(chromosome_names[[i]]))), "...")

    plus_matches  <- matchPattern(DS_sequence, chromosome_sequences_list[[i]], max.mismatch = 1)
    minus_matches <- matchPattern(DS_sequence_reversed, chromosome_sequences_list[[i]], max.mismatch = 1)
    if (length(plus_matches) > 0) {
      results_list <- c(results_list, list(XStringViewsToDf(plus_matches, chromosome_names[[i]], "+", my_sequence)))
    }
    if (length(minus_matches) > 0) {
      results_list <- c(results_list, list(XStringViewsToDf(minus_matches, chromosome_names[[i]], "-", my_sequence)))
    }
  }

  if (length(results_list) == 0) {
    return(NULL)
  } else {
    results_df <- do.call(rbind.data.frame, c(results_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
    return(results_df)
  }
}


HarmonizeSequenceLengths <- function(sequences_vec, verbose = TRUE) {
  string_lengths <- vapply(sequences_vec, nchar, integer(1))
  use_length <- min(string_lengths)
  are_longer <- string_lengths > use_length
  if (any(are_longer)) {
    truncation_string <- paste0(" ", sum(are_longer), " sequence", if (sum(are_longer) > 1) "s" else "",
                                " were longer than ", use_length, " bp and were truncated."
                                )
    sequences_vec <- substr(sequences_vec, 1, use_length)
  } else {
    truncation_string <- ""
  }
  if (verbose) {
    message(paste0(use_length, "-mer sequences will be used.", truncation_string))
  }
  results_list <- list("Truncated_sequences" = sequences_vec, "Length" = use_length)
  return(results_list)
}



RetrieveSequence <- function(chromosome, strand, start, end) {
  result_string <- as.character(Views(chromosome_sequences_list[[chromosome]], start, end))
  if (strand == "-") {
    result_string <- as.character(reverseComplement(DNAString(result_string)))
  } else if (strand != "+") {
    stop("Invalid value for the strand parameter!")
  }
  return(result_string)
}



FindVariableLengthSequences <- function(sequences_vec) {
  results_list <- lapply(sequences_vec, FindSequence)
  results_df <- do.call(rbind.data.frame, c(results_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  return(results_df)
}



FindSequences <- function(sequences_vec) {

  num_sequences <- length(sequences_vec)
  message(paste0("Looking for matches for ", num_sequences, " sequence", if (num_sequences > 1) "s" else "", " in the Hsapiens.UCSC.hg38 genome..."))

  harmonized_sequences <- HarmonizeSequenceLengths(sequences_vec)
  sequences_vec <- toupper(harmonized_sequences[["Truncated_sequences"]])

  my_PDict          <- PDict(sequences_vec, max.mismatch = 1)
  my_PDict_reversed <- PDict(reverseComplement(DNAStringSet(sequences_vec)), max.mismatch = 1)

  results_list <- vector("list", length(chromosome_names))
  names(results_list) <- chromosome_names

  for (i in seq_along(chromosome_names)) {

    message(paste0("Searching the + strand of chromosome ", substr(chromosome_names[[i]], 4, nchar(chromosome_names[[i]]))), "...")
    plus_ranges_list <- matchPDict(my_PDict, chromosome_sequences_list[[i]], max.mismatch = 1)

    message(paste0("Searching the - strand of chromosome ", substr(chromosome_names[[i]], 4, nchar(chromosome_names[[i]]))), "...")
    minus_ranges_list <- matchPDict(my_PDict_reversed, chromosome_sequences_list[[i]], max.mismatch = 1)

    plus_matches_df_list <- lapply(seq_along(sequences_vec), function(x) {
      if (length(plus_ranges_list[[x]]) == 0) {
        return(NULL)
      } else {
        return(XStringViewsToDf(Views(chromosome_sequences_list[[i]], start(plus_ranges_list[[x]]), end(plus_ranges_list[[x]])),
                                chromosome_names[[i]], "+", sequences_vec[[x]]
                                )
               )
      }
    })
    minus_matches_df_list <- lapply(seq_along(sequences_vec), function(x) {
      if (length(minus_ranges_list[[x]]) == 0) {
        return(NULL)
      } else {
        return(XStringViewsToDf(Views(chromosome_sequences_list[[i]], start(minus_ranges_list[[x]]), end(minus_ranges_list[[x]])),
                                chromosome_names[[i]], "-", sequences_vec[[x]]
                                )
               )
      }
    })

    results_df <- do.call(rbind.data.frame, c(plus_matches_df_list, minus_matches_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
    if (nrow(results_df) > 0) {
      results_list[[i]] <- results_df
    } else {
      results_list[[i]] <- NULL
    }
  }

  combined_results_df <- do.call(rbind.data.frame, c(results_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))

  new_order <- order(match(combined_results_df[["Reference"]], sequences_vec),
                     match(combined_results_df[["Chromosome"]], chromosome_names),
                     combined_results_df[["Start"]],
                     match(combined_results_df[["Strand"]], c("+", "-")),
                     combined_results_df[["End"]]
                     )

  combined_results_df <- combined_results_df[new_order, ]
  row.names(combined_results_df) <- NULL

  return(combined_results_df)
}




