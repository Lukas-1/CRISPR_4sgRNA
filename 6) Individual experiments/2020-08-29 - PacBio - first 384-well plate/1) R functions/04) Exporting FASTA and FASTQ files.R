### 20th October 2020 ###


# Import packages and source code -----------------------------------------

library("Biostrings")



# Define functions --------------------------------------------------------

FormatFixedWidthInteger <- function(integer_vec, full_vec = integer_vec) {
  integer_width <- max(nchar(as.character(as.integer(full_vec))))
  result <- formatC(integer_vec, width = integer_width, flag = "0")
  return(result)
}



BuildChunksDf <- function(total_number, number_per_file = 20L) {
  entry_vec <- seq_len(total_number)
  num_files <- ceiling(total_number / number_per_file)
  file_vec <- rep(seq_len(num_files), each = number_per_file)
  file_vec <- file_vec[seq_len(total_number)]

  first_file <- tapply(entry_vec, file_vec, min)
  last_file <- tapply(entry_vec, file_vec, max)

  chunk_range <- paste0(FormatFixedWidthInteger(first_file, full_vec = last_file),
                        "_to_",
                        FormatFixedWidthInteger(last_file)
                        )
  results_df <- data.frame("Entry_number"   = seq_len(total_number),
                           "File_number"    = file_vec,
                           "File_name"      = chunk_range[file_vec],
                           stringsAsFactors = FALSE,
                           row.names        = NULL
                           )
  return(results_df)
}



NoteLongReads <- function(reads_object, are_too_long, export_dir, file_name, max_length) {
  assign("delete_reads_object", reads_object, envir = globalenv())
  assign("delete_are_too_long", are_too_long, envir = globalenv())
  stopifnot(!(all(are_too_long)))
  if (any(are_too_long)) {
    long_IDs <- names(reads_object)[are_too_long]
    write_message <- paste0("The following reads were too long (>",
                            max_length, " bp) and could not",
                            " be exported: ", paste0(long_IDs, collapse = ", ")
                            )
    message(write_message)
    write.table(write_message,
                file = file.path(export_dir, paste0(file_name, " - log.txt")),
                quote = FALSE, row.names = FALSE, col.names = FALSE
                )
  }
}



ExportFASTA <- function(export_fasta, export_dir, file_name, max_length = 20000L) {
  are_too_long <- lengths(export_fasta) > max_length
  NoteLongReads(export_fasta, are_too_long, export_dir, file_name, max_length)
  writeXStringSet(export_fasta[!(are_too_long)],
                  filepath = file.path(export_dir, paste0(file_name, ".fasta"))
                  )
  return(invisible(NULL))
}


ExportFASTQ <- function(export_fastq, export_dir, file_name, max_length = 20000L) {
  are_too_long <- lengths(export_fastq) > max_length
  NoteLongReads(export_fastq, are_too_long, export_dir, file_name, max_length)
  writeQualityScaledXStringSet(export_fastq[!(are_too_long)],
                               filepath = file.path(export_dir, paste0(file_name, ".fastq"))
                               )
  return(invisible(NULL))
}



MakeCategoryString <- function(ccs_df, extracted_df) {
  abbreviations_vec <- c(
    "Correct"            = "P",
    "Flanking insertion" = "P",
    "Mutation"           = "M",
    "Contamination"      = "C",
    "Deletion"           = "D"
  )
  are_present <- ccs_df[["ZMW"]] %in% extracted_df[["ZMW"]]
  features_mat <- GetFeaturesData(extracted_df, ccs_df[["ZMW"]][are_present])
  cat_mat <- features_mat[, paste0("sg", 1:4, "_cr", 1:4, "_category")]
  for (i in 1:4) {
    cat_mat[, i] <- abbreviations_vec[cat_mat[, i]]
  }
  cat_vec <- do.call(paste, c(lapply(1:4, function(i) paste0(cat_mat[, i], i)), sep = "_"))
  cat_vec <- ifelse(cat_vec == "P1_P2_P3_P4", "4_correct", cat_vec)
  results_vec <- rep(NA, nrow(ccs_df))
  results_vec[are_present] <- cat_vec
  return(results_vec)
}



ExportSequences <- function(ccs_df,
                            fasta_output_dir,
                            fastq_output_dir,
                            append_to_file_name = "",
                            prefer_ccs          = TRUE,
                            use_zmws            = NULL,
                            split_into_chunks   = FALSE,
                            chunk_size          = 50L,
                            ID_column           = "Well_number",
                            unique_IDs          = seq_len(384),
                            export_fasta        = TRUE
                            ) {

  # Ignore wells with zero reads
  ID_matches <- match(unique_IDs, ccs_df[[ID_column]])
  are_present <- !(is.na(ID_matches))
  unique_IDs <- unique_IDs[are_present]
  ID_matches <- ID_matches[are_present]

  if (ID_column == "Well_number") {
    wells_formatted <- formatC(unique_IDs, flag = "0", width = 3)
    well_names <- paste0("well", wells_formatted)
  } else {
    wells_formatted <- unique_IDs
    well_names <- unique_IDs
    has_plate_number <- "Plate_number" %in% names(ccs_df)
    if (has_plate_number) {
      plate_numbers <- ccs_df[["Plate_number"]][ID_matches]
      plate_names <- paste0("Plate", formatC(plate_numbers, width = 2, flag = "0"))
    }
  }

  file_names <- paste0(well_names, append_to_file_name)

  lima_zmws <- ccs_df[["ZMW"]][ccs_df[["Passed_filters"]]]
  if (!(is.null(use_zmws))) {
    stopifnot(all(use_zmws %in% lima_zmws))
    lima_zmws <- use_zmws
  }

  ccs_matches_vec <- match(lima_zmws, ccs_df[["ZMW"]])

  lima_well_IDs <- ccs_df[ccs_matches_vec, ID_column]
  ZMW_strings <- ccs_df[ccs_matches_vec, "ZMW_string"]

  for (i in seq_along(unique_IDs)) {

    are_this_well <- lima_well_IDs %in% unique_IDs[[i]]
    any_reads <- any(are_this_well)

    if (any_reads) {
      this_well_zmws <- lima_zmws[are_this_well]
      ccs_matches <- match(this_well_zmws, ccs_df[["ZMW"]])
      stopifnot(!(anyNA(ccs_matches)))
      if (prefer_ccs) {
        export_seq <- ccs_df[["Sequence"]][ccs_matches]
        export_qual <- ccs_df[["Quality"]][ccs_matches]
      } else {
        export_seq <- substr(ccs_df[["Sequence"]][ccs_matches],
                             ccs_df[["Clip_start"]][ccs_matches],
                             ccs_df[["Clip_end"]][ccs_matches]
                             )
        export_qual <- substr(ccs_df[["Quality"]][ccs_matches],
                              ccs_df[["Clip_start"]][ccs_matches],
                              ccs_df[["Clip_end"]][ccs_matches]
                              )
      }
      names(export_seq) <- ZMW_strings[are_this_well]
      export_seq <- DNAStringSet(export_seq)
      export_fastq <- QualityScaledBStringSet(export_seq, PhredQuality(export_qual))
      message(paste0("Exporting reads for well ", wells_formatted[[i]], "..."))
    } else {
      message(paste0("No reads were present for well ", wells_formatted[[i]], "!"))
    }

    if ("Plate_number" %in% names(ccs_df)) {
      plate_folder <- plate_names[[i]]
      fasta_dir <- file.path(fasta_output_dir, plate_folder)
      fastq_dir <- file.path(fastq_output_dir, plate_folder)
      if (export_fasta) {
        dir.create(fasta_dir, showWarnings = FALSE)
      }
      dir.create(fastq_dir, showWarnings = FALSE)
    } else {
      fasta_dir <- fasta_output_dir
      fastq_dir <- fastq_output_dir
    }

    if (split_into_chunks) {
      fasta_folder <- file.path(fasta_dir, well_names[[i]])
      fastq_folder <- file.path(fastq_dir, well_names[[i]])
      if (export_fasta) {
        dir.create(fasta_folder, showWarnings = FALSE)
      }
      dir.create(fastq_folder, showWarnings = FALSE)
      if (!(any_reads)) {
        next
      }
      num_reads <- length(this_well_zmws)
      chunks_df <- BuildChunksDf(num_reads, chunk_size)
      for (file_number in unique(chunks_df[["File_number"]])) {
        are_this_file <- chunks_df[["File_number"]] == file_number
        file_name <- paste0(file_names[[i]], "__", unique(chunks_df[["File_name"]][are_this_file]))
        if (export_fasta) {
          ExportFASTA(export_seq[are_this_file], fasta_folder, file_name)
        }
        ExportFASTQ(export_fastq[are_this_file], fastq_folder, file_name)
      }
    } else {
      if (!(any_reads)) {
        next
      }
      if (export_fasta) {
        ExportFASTA(export_seq, fasta_dir, file_names[[i]])
      }
      ExportFASTQ(export_fastq, fastq_dir, file_names[[i]])
    }
  }
  if (split_into_chunks) {
    message("")
  }
}


