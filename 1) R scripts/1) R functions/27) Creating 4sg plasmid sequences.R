### 27th June 2020 ###





# Import packages and source code -----------------------------------------

library("Biostrings")




# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
plasmids_directory <- file.path(CRISPR_root_directory, "2) Input data", "Vectors")





# Read in data ------------------------------------------------------------

plasmid_genbank_snapgene <- scan(file = file.path(plasmids_directory, "pYJA5 empty vector for 4sg expression_04032019.gbk"),
                                 what = character(),
                                 sep = "\n"
                                 )
first_index <- which(plasmid_genbank_snapgene == "ORIGIN") + 1L
prior_lines <- plasmid_genbank_snapgene[seq_len(first_index - 1)]
last_index <- first_index + 159
after_seq <- seq(from = last_index + 1,
                 to = length(plasmid_genbank_snapgene)
                 )
after_lines <- plasmid_genbank_snapgene[after_seq]
sequence_seq <- seq(from = first_index, to = last_index)
sequence_lines <- plasmid_genbank_snapgene[sequence_seq]





# Define functions --------------------------------------------------------

ReplaceNNNLines <- function(lines_vec, four_guides) {

  ### Extract the positions of the 4 guides

  matches_list <- gregexpr("n+[n ]*", lines_vec)
  are_matches <- matches_list != -1
  lines_vec[are_matches]


  ### Group locations that span two lines

  guide_numbers <- rep(NA_integer_, sum(are_matches))
  current_number <- 1L
  current_char_sum <- 0L
  for (i in seq_along(guide_numbers)) {
    this_index <- which(are_matches)[[i]]
    guide_numbers[[i]] <- current_number
    current_char_sum <- current_char_sum + attributes(matches_list[[which(are_matches)[[i]]]])[["match.length"]]
    if (current_char_sum >= 21) {
      current_char_sum <- 0
      current_number <- current_number + 1L
    } else {
      stopifnot((which(are_matches)[[i + 1]] - this_index) == 1)
    }
  }

  ### Prepare for replacements

  are_4sg <- guide_numbers %in% 1:4  # The 5th location seems to be something else!
  are_matches[are_matches] <- are_4sg
  guide_numbers <- guide_numbers[are_4sg]

  only_matches_list_list <- split(matches_list[are_matches], guide_numbers)
  indices_list <- split(seq_len(sum(are_matches)), guide_numbers)

  new_lines_vec <- rep(NA_character_, sum(are_matches))
  sequences_vec <- rev(four_guides)
  sequences_vec <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(sequences_vec)))


  ### Replace the 4 guide sequences

  for (i in seq_along(only_matches_list_list)) {
    sum_replaced <- 0L
    for (j in indices_list[[i]]) {
      match_index <- which(are_matches)[[j]]
      nnn_start <- matches_list[[match_index]][[1]]
      nnn_length <- attributes(matches_list[[match_index]])[["match.length"]]
      nnn_end <- nnn_start + nnn_length - 1L
      nnn_string <- substr(lines_vec[[match_index]], nnn_start, nnn_end)
      nnn_chars <- strsplit(nnn_string, "", fixed = TRUE)[[1]]
      replacement_chars <- rep(NA_character_, nnn_length)
      for (char_index in seq_len(nnn_length)) {
        if (nnn_chars[[char_index]] == " ") {
          replacement_chars[[char_index]] <- " "
        } else {
          sum_replaced <- sum_replaced + 1L
          assign("delete_sequences_vec", sequences_vec, envir = globalenv())
          replacement_chars[[char_index]] <- substr(sequences_vec[[i]],
                                                    sum_replaced,
                                                    sum_replaced
                                                    )
        }
      }
      replacement_string <- paste0(replacement_chars, collapse = "")
      if (nnn_start != 1) {
        string_before <- substr(lines_vec[[match_index]], 1, nnn_start - 1L)
      } else {
        string_before <- ""
      }
      line_length <- nchar(lines_vec[[match_index]])
      if (nnn_end != line_length) {
        string_after <- substr(lines_vec[[match_index]], nnn_end + 1L, line_length)
      } else {
        string_after <- ""
      }
      new_lines_vec[[j]] <- paste0(string_before, replacement_string, string_after)
    }
    stopifnot(sum_replaced == 20)
  }

  results_vec <- lines_vec
  results_vec[are_matches] <- new_lines_vec
  return(results_vec)
}




ExportVectorsForGene <- function(symbol_or_entrez, CRISPR_df) {

  if (!(is.na(suppressWarnings(as.integer(symbol_or_entrez))))) {
    are_this_gene <- CRISPR_df[["Entrez_ID"]] %in% as.character(symbol_or_entrez)
  } else {
    are_this_gene <- CRISPR_df[["Gene_symbol"]] %in% symbol_or_entrez
  }
  are_top4 <- CRISPR_df[["Rank"]] %in% 1:4
  are_chosen <- are_this_gene & are_top4
  num_guides <- sum(are_chosen)
  num_TSSs <- num_guides / 4

  if (num_guides == 0) {
    stop(paste0("The gene '", symbol_or_entrez, "' was not found!"))
  }
  stopifnot((num_guides %% 4) == 0)

  if (num_TSSs > 1) {
    file_name_postfixes <- paste0("_TSS", seq_len(num_TSSs))
  } else {
    file_name_postfixes <- ""
  }
  sequences_list <- split(CRISPR_df[["sgRNA_sequence"]][are_chosen],
                          rep(seq_len(num_TSSs), each = 4)
                          )
  for (i in seq_along(sequences_list)) {
    new_lines <- ReplaceNNNLines(sequence_lines, sequences_list[[i]])
    file_name <- paste0("4sg_inserted_vector_",
                        symbol_or_entrez,
                        file_name_postfixes[[i]],
                        ".gbk"
                        )
    write.table(x         = c(prior_lines, new_lines, after_lines),
                file      = file.path(plasmid_output_directory, file_name),
                col.names = FALSE,
                row.names = FALSE,
                quote     = FALSE
                )
  }
  return(invisible(NULL))
}








