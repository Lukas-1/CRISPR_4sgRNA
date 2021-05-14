### 28th October 2020 ###




# Import packages and source code -----------------------------------------

library("Biostrings")





# Define functions --------------------------------------------------------

ReadsForZMW <- function(zmw, ccs_df, use_fastq = TRUE) {

  stopifnot(all(c("subreads_bam", "subreads_stats_df") %in% ls(envir = globalenv())))

  subread_zmws <- subreads_stats_df[["hole"]]

  are_ccs_zmw <- ccs_df[["ZMW"]] == zmw
  stopifnot(sum(are_ccs_zmw) == 1)
  are_subread_zmw <- subread_zmws == zmw

  all_sequences <- c(DNAStringSet(ccs_df[["Sequence"]][are_ccs_zmw]),
                     subreads_bam[["seq"]][are_subread_zmw]
                     )
  starts_vec <- subreads_stats_df[["start"]][are_subread_zmw]
  ends_vec <- subreads_stats_df[["end"]][are_subread_zmw]
  names(all_sequences) <- paste0(zmw, "_", c("ccs", paste0(starts_vec, "_", ends_vec)))

  if (use_fastq) {
    all_qualities <- c(PhredQuality(ccs_df[["Quality"]][are_ccs_zmw]),
                       rep("!", sum(are_subread_zmw))
                       )
    results_object <- QualityScaledBStringSet(all_sequences, all_qualities)
  } else {
    results_object <- all_sequences
  }
  return(results_object)
}





ExportReadsForZMW <- function(zmw,
                              ccs_df,
                              output_dir = subreads_zmws_directory,
                              use_fastq = TRUE,
                              subread_max_length = 20000L
                              ) {
  are_ccs_zmw <- ccs_df[["ZMW"]] == zmw

  stopifnot(sum(are_ccs_zmw) == 1)
  zmw_well_number <- ccs_df[["Well_number"]][are_ccs_zmw]

  if (is.na(zmw_well_number)) {
    well_prefix <- "no_well"
  } else {
    well_prefix <- paste0("well", formatC(zmw_well_number, flag = "0", width = 3))
  }
  file_name <- paste0(well_prefix, "_", zmw)

  reads_object <- ReadsForZMW(zmw, ccs_df, use_fastq = use_fastq)
  are_too_long <- lengths(reads_object) > subread_max_length
  stopifnot(!(all(are_too_long)))
  if (any(are_too_long)) {
    long_IDs <- names(reads_object)[are_too_long]
    write_message <- paste0("The following reads were too long (>",
                            subread_max_length, " bp) and could not",
                            " be exported: ", paste0(long_IDs, collapse = ", ")
                            )
    if (!(use_fastq)) {
      message(write_message)
    }
    write.table(write_message,
                file = file.path(output_dir, paste0(file_name, " - log.txt")),
                quote = FALSE, row.names = FALSE, col.names = FALSE
                )
  }
  if (use_fastq) {
    writeQualityScaledXStringSet(reads_object[!(are_too_long)],
                                 filepath = file.path(output_dir, paste0(file_name, ".fastq"))
                                 )
  } else {
    writeXStringSet(reads_object[!(are_too_long)],
                    filepath = file.path(output_dir, paste0(file_name, ".fasta"))
                    )
  }

  return(invisible(NULL))
}





ExportSubreadsForWells <- function(ccs_df,
                                   fasta_output_dir,
                                   fastq_output_dir,
                                   prefer_ccs   = TRUE,
                                   use_zmws     = NULL,
                                   max_num_ZMWs = 20L,
                                   wells_vec    = seq_len(384)
                                   ) {

  wells_formatted <- formatC(seq_len(384), flag = "0", width = 3)
  well_names <- paste0("well", wells_formatted)

  for (i in wells_vec) {
    message(paste0("Exporting reads for well ", wells_formatted[[i]], "..."))
    are_this_well <- (ccs_df[["Well_number"]] %in% i) &
                     (ccs_df[["Passed_filters"]] == 1)
    this_well_zmws <- ccs_df[["ZMW"]][are_this_well]
    fasta_folder <- file.path(fasta_output_dir, well_names[[i]])
    fastq_folder <- file.path(fastq_output_dir, well_names[[i]])
    dir.create(fasta_folder, showWarnings = FALSE)
    dir.create(fastq_folder, showWarnings = FALSE)
    if (length(this_well_zmws) > max_num_ZMWs) {
      this_well_zmws <- this_well_zmws[seq_len(max_num_ZMWs)]
    }
    for (zmw in this_well_zmws) {
      # message(paste0("Exporting reads for ZMW ", zmw))
      ExportReadsForZMW(zmw, ccs_df, output_dir = fasta_folder, use_fastq = FALSE)
      ExportReadsForZMW(zmw, ccs_df, output_dir = fastq_folder, use_fastq = TRUE)
    }
  }
  message("")
}



