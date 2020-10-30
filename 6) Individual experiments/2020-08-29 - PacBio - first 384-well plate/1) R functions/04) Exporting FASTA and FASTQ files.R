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



ExportSequences <- function(lima_reads,
                            report_df,
                            fasta_output_dir,
                            fastq_output_dir,
                            ccs_reads           = NULL,
                            append_to_file_name = "",
                            prefer_ccs          = TRUE,
                            use_zmws            = NULL,
                            split_into_chunks   = FALSE,
                            chunk_size          = 50L,
                            wells_vec           = seq_len(384)
                            ) {

  wells_formatted <- formatC(seq_len(384), flag = "0", width = 3)
  well_names <- paste0("well", wells_formatted)
  file_names <- paste0(well_names, append_to_file_name)

  lima_zmws <- as.integer(substr(lima_reads[["qname"]], 22, nchar(lima_reads[["qname"]]) - 4))
  if (!(is.null(use_zmws))) {
    stopifnot(all(use_zmws %in% lima_zmws))
    lima_zmws <- use_zmws
  }

  ccs_zmws <- as.integer(substr(report_df[["ZMW"]], 22, nchar(report_df[["ZMW"]])))
  ccs_well_numbers <- GetWellNumbers(report_df)
  lima_well_numbers <- ccs_well_numbers[match(lima_zmws, ccs_zmws)]
  use_ccs <- prefer_ccs && (!(is.null(ccs_reads)))
  if (!(is.null(ccs_reads))) {
    stopifnot(identical(ccs_zmws, as.integer(substr(ccs_reads[["qname"]], 22, nchar(ccs_reads[["qname"]]) - 4))))
  }

  for (i in wells_vec) {
    are_this_well <- lima_well_numbers %in% i
    this_well_zmws <- lima_zmws[are_this_well]
    if (prefer_ccs && !(is.null(ccs_reads))) {
      ccs_matches <- match(this_well_zmws, ccs_zmws)
      stopifnot(!(anyNA(ccs_matches)))
      export_seq <- ccs_reads[["seq"]][ccs_matches]
      export_qual <- ccs_reads[["qual"]][ccs_matches]
    } else {
      export_seq <- lima_reads[["seq"]][are_this_well]
      export_qual <- lima_reads[["qual"]][are_this_well]
    }
    names(export_seq) <- as.character(this_well_zmws)
    export_fastq <- QualityScaledBStringSet(export_seq, export_qual)

    if (split_into_chunks) {
      message(paste0("Exporting reads for well ", wells_formatted[[i]]), "...")
      num_reads <- length(this_well_zmws)
      chunks_df <- BuildChunksDf(num_reads, chunk_size)
      fasta_folder <- file.path(fasta_output_dir, well_names[[i]])
      fastq_folder <- file.path(fastq_output_dir, well_names[[i]])
      dir.create(fasta_folder, showWarnings = FALSE)
      dir.create(fastq_folder, showWarnings = FALSE)
      for (file_number in chunks_df[["File_number"]]) {
        are_this_file <- chunks_df[["File_number"]] == file_number
        file_name <- paste0(file_names[[i]], "_", unique(chunks_df[["File_name"]][are_this_file]))
        writeXStringSet(export_seq[are_this_file],
                        filepath = file.path(fasta_folder, paste0(file_name, ".fasta"))
                        )
        writeQualityScaledXStringSet(export_fastq[are_this_file],
                                     filepath = file.path(fastq_folder, paste0(file_name, ".fastq"))
                                     )
      }
    } else {
      fasta_path <- file.path(fasta_output_dir, paste0(file_names[[i]], ".fasta"))
      fastq_path <- file.path(fastq_output_dir, paste0(file_names[[i]], ".fastq"))
      writeXStringSet(export_seq, filepath = fasta_path)
      writeQualityScaledXStringSet(export_fastq, filepath = fastq_path)
    }
  }
  if (split_into_chunks) {
    message("")
  }
}


