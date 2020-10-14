### 9th October 2020 ###





# Import packages and source code -----------------------------------------

library("Biostrings")




# Define folder paths -----------------------------------------------------

CRISPR_root_directory       <- "~/CRISPR"
file_directory              <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
file_input_directory        <- file.path(file_directory, "2) Input")
R_objects_directory         <- file.path(file_directory, "3) R objects")

metadata_directory          <- file.path(file_input_directory, "Metadata")
sg_sequences_file           <- file.path(metadata_directory, "1-384 sgRNA summary.xlsm")
tracRNAs_file               <- file.path(metadata_directory, "tracRNAs_4sg.txt")
promoters_file              <- file.path(metadata_directory, "promoters_4sg.txt")
full_plasmid_file           <- file.path(metadata_directory, "reference with 4sg and barcode example.fa")

annotated_plasmid_directory <- file.path(file_input_directory, "Annotated plasmids")

file_output_directory       <- file.path(file_directory, "5) Output")
reference_output_directory  <- file.path(file_output_directory, "Reference sequences")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "1) Process and export barcodes.RData"))
load(file.path(R_objects_directory, "3) Import and process sgRNA sequences.RData"))





# Read in data ------------------------------------------------------------

tracRNAs_vec <- read.table(tracRNAs_file, header = FALSE,
                           stringsAsFactors = FALSE
                           )[[1]]
promoters_vec <- read.table(promoters_file, header = FALSE,
                            stringsAsFactors = FALSE
                            )[[1]]
plasmid_lines_vec <- read.table(full_plasmid_file, header = FALSE, skip = 1,
                                stringsAsFactors = FALSE
                                )[[1]]




# Define functions --------------------------------------------------------


ReplaceNNNLines <- function(lines_vec, sequences_vec) {

  assign("delete_lines_vec", lines_vec, envir = globalenv())
  assign("delete_sequences_vec", sequences_vec, envir = globalenv())

  ### Extract the positions of the 4 guides

  matches_list <- gregexpr("n+[n ]*", lines_vec)
  are_matches <- matches_list != -1


  ### Group locations that span two lines

  match_numbers <- rep(NA_integer_, sum(are_matches))
  current_number <- 0L
  for (i in seq_along(match_numbers)) {
    if (i == 1) {
      previous_reaches_line_end <- FALSE
    } else {
      previous_reaches_line_end <- grepl("n+[n ]*$", lines_vec[are_matches][[i - 1]])
    }
    if (previous_reaches_line_end) {
      starts_from_line_start <- grepl("^[n 0-9]*n", lines_vec[are_matches][[i]])
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

  only_matches_list_list <- split(matches_list[are_matches], match_numbers)
  indices_list <- split(seq_len(sum(are_matches)), match_numbers)

  new_lines_vec <- rep(NA_character_, sum(are_matches))


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



ReadInPlasmids <- function(file_path) {

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


PlasmidForWell <- function(well_number) {
  stopifnot("plasmid_list" %in% ls(envir = globalenv()))
  search_string <- paste0("_well", well_number, ".gbk")
  plasmid_name <- grep(search_string, names(plasmid_list), fixed = TRUE)
  if (length(plasmid_name) == 0) {
    plasmid_name <- "pYJA5_with_4sg_20bp.gbk"
  }
  return(plasmid_name)
}





# Read in annotated plasmid sequences -------------------------------------

plasmid_files <- list.files(annotated_plasmid_directory)

plasmid_list <- sapply(plasmid_files, function(x) {
  ReadInPlasmids(file.path(annotated_plasmid_directory, x))
}, simplify = FALSE)






# Define the reference sequences ------------------------------------------

rev_tracRNAs_vec <- as.character(reverseComplement(DNAStringSet(tracRNAs_vec)))

guides_ref_list <- lapply(1:4, function(x) {
  sg_column <- paste0("Sequence_sg", x)
  paste0(sg_sequences_df[[sg_column]], rev_tracRNAs_vec[[x]], "TTTTT")
})

guides_with_promoters_list <- lapply(1:4, function(x) {
  paste0(toupper(promoters_vec[[x]]), guides_ref_list[[x]])
})

plasmid_string <- toupper(paste0(plasmid_lines_vec, collapse = ""))
plasmid_string <- substr(plasmid_string, 11, nchar(plasmid_string) - 10)
plasmids_vec <- vapply(seq_len(384), function(x) {
  sg_N <- paste0(rep("N", 20), collapse = "")
  sg_sequences <- vapply(1:4, function(y) {
    sg_sequences_df[[paste0("Sequence_sg", y)]][[x]]
  }, "")
  result_string <- sub(sg_N, sg_sequences[[1]], plasmid_string, fixed = TRUE)
  result_string <- sub(sg_N, sg_sequences[[2]], result_string, fixed = TRUE)
  result_string <- sub(sg_N, sg_sequences[[3]], result_string, fixed = TRUE)
  result_string <- sub(sg_N, sg_sequences[[4]], result_string, fixed = TRUE)
  return(result_string)
}, "")






# Insert barcodes and guides into the annotated plasmids ------------------

plasmid_lines_list <- lapply(seq_len(384), function(x) {
  use_plasmid <- PlasmidForWell(x)
  use_lines <- plasmid_list[[use_plasmid]][["sequence"]]
  sg_seqs <- as.character(sg_sequences_df[x, paste0("Sequence_sg", 1:4)])
  replace_seqs <- c(column_bc_vec[[x]], sg_seqs, row_bc_vec[[x]])
  ReplaceNNNLines(use_lines, replace_seqs)
})





# Export the annotated plasmids -------------------------------------------

for (i in seq_len(384)) {
  file_name <- paste0("Well",
                      formatC(i, flag = "0", width = 3),
                      "_barcoded_cassette.gbk"
                      )
  use_plasmid <- PlasmidForWell(i)
  export_lines <- c(plasmid_list[[use_plasmid]][["preceding"]],
                    plasmid_lines_list[[i]],
                    plasmid_list[[use_plasmid]][["following"]]
                    )
  write.table(x         = export_lines,
              file      = file.path(reference_output_directory, file_name),
              col.names = FALSE,
              row.names = FALSE,
              quote     = FALSE
              )
}






# Map barcode IDs to wells ------------------------------------------------

all_barcodes_df <- expand.grid(paste0("RevBC", seq_len(24)),
                               paste0("FwdBC", seq_len(16)),
                               KEEP.OUT.ATTRS = FALSE,
                               stringsAsFactors = FALSE
                               )[, 2:1]
barcodes_to_wells_map <- seq_len(384)
names(barcodes_to_wells_map) <- paste0(all_barcodes_df[[1]],
                                       "--",
                                       all_barcodes_df[[2]]
                                       )






# Save data ---------------------------------------------------------------

save(list = c("guides_ref_list", "guides_with_promoters_list",
              "plasmid_string", "plasmids_vec"
              ),
     file = file.path(R_objects_directory, "4) Create reference sequences for each well - raw sequences.RData")
     )

save(list = "plasmid_list",
     file = file.path(R_objects_directory, "4) Create reference sequences for each well - annotated plasmid.RData")
     )









