### 28th October 2020 ###




# Import packages and source code -----------------------------------------

library("Biostrings")





# Define functions --------------------------------------------------------

ReadsForZMW <- function(zmw, use_sl7 = TRUE, use_fastq = use_fastq) {

  stopifnot(all(c("subreads_bam", "subreads_stats_df") %in% ls(envir = globalenv())))

  if (use_sl7) {
    ccs_list <- sl7_ccs3_ccs
  } else {
    ccs_list <- sl9_ccs3_ccs
  }

  ccs_zmws <- as.integer(substr(ccs_list[["qname"]], 22, nchar(ccs_list[["qname"]]) - 4))
  subread_zmws <- subreads_stats_df[["hole"]]

  are_ccs_zmw <- ccs_zmws == zmw
  stopifnot(sum(are_ccs_zmw) == 1)
  are_subread_zmw <- subread_zmws == zmw

  all_sequences <- c(ccs_list[["seq"]][are_ccs_zmw],
                     subreads_bam[["seq"]][are_subread_zmw]
                     )
  starts_vec <- subreads_stats_df[["start"]][are_subread_zmw]
  ends_vec <- subreads_stats_df[["end"]][are_subread_zmw]
  names(all_sequences) <- paste0(zmw, "_", c("ccs", paste0(starts_vec, "_", ends_vec)))

  if (use_fastq) {
    all_qualities <- c(ccs_list[["qual"]][are_ccs_zmw],
                       rep("!", sum(are_subread_zmw))
                       )
    results_object <- QualityScaledBStringSet(all_sequences, all_qualities)
  } else {
    results_object <- all_sequences
  }
  return(results_object)
}





ExportReadsForZMW <- function(zmw,
                              use_sl7 = TRUE,
                              output_dir = subreads_zmws_directory,
                              use_fastq = TRUE
                              ) {
  if (use_sl7) {
    report_df <- sl7_ccs3_report_df
  } else {
    report_df <- sl9_ccs3_report_df
  }
  ccs_well_numbers <- GetWellNumbers(report_df)
  ccs_zmws <- as.integer(substr(report_df[["ZMW"]], 22, nchar(report_df[["ZMW"]])))
  are_ccs_zmw <- ccs_zmws == zmw

  stopifnot(sum(are_ccs_zmw) == 1)
  zmw_well_number <- ccs_well_numbers[are_ccs_zmw]

  if (is.na(zmw_well_number)) {
    well_prefix <- "no_well"
  } else {
    well_prefix <- paste0("well", formatC(zmw_well_number, flag = "0", width = 3))
  }
  file_name <- paste0(well_prefix, "_", zmw)

  reads_object <- ReadsForZMW(zmw, use_sl7 = use_sl7, use_fastq = use_fastq)
  are_too_long <- lengths(reads_object) > 25000
  stopifnot(!(all(are_too_long)))
  if (any(are_too_long)) {
    if (any(are_too_long)) {
      long_IDs <- names(reads_object)[are_too_long]
      write_message <- paste0("The following reads were too long and could not",
                              " be exported: ", paste0(long_IDs, collapse = ", ")
                              )
      write.table(write_message,
                  file = file.path(output_dir, paste0(file_name, " - log.txt")),
                  quote = FALSE, row.names = FALSE, col.names = FALSE
                  )
    }
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





ExportSubreadsForWells <- function(fasta_output_dir,
                                   fastq_output_dir,
                                   prefer_ccs   = TRUE,
                                   use_zmws     = NULL,
                                   max_num_ZMWs = 20L,
                                   wells_vec    = seq_len(384),
                                   use_sl7      = TRUE
                                   ) {

  wells_formatted <- formatC(seq_len(384), flag = "0", width = 3)
  well_names <- paste0("well", wells_formatted)

  if (use_sl7) {
    lima_list <- sl7_ccs3_lima
    ccs_list <- sl7_ccs3_ccs
    report_df <- sl7_ccs3_report_df
  } else {
    lima_list <- sl7_ccs3_lima
    ccs_list <- sl7_ccs3_ccs
    report_df <- sl9_ccs3_report_df
  }

  lima_zmws <- as.integer(substr(lima_list[["qname"]], 22, nchar(lima_list[["qname"]]) - 4))

  if (!(is.null(use_zmws))) {
    stopifnot(all(use_zmws %in% lima_zmws))
    lima_zmws <- use_zmws
  }

  ccs_zmws <- as.integer(substr(report_df[["ZMW"]], 22, nchar(report_df[["ZMW"]])))
  ccs_well_numbers <- GetWellNumbers(report_df)
  lima_well_numbers <- ccs_well_numbers[match(lima_zmws, ccs_zmws)]

  stopifnot(identical(ccs_zmws, as.integer(substr(ccs_list[["qname"]], 22, nchar(ccs_list[["qname"]]) - 4))))

  for (i in wells_vec) {
    message(paste0("Exporting reads for well ", wells_formatted[[i]], "..."))
    are_this_well <- lima_well_numbers %in% i
    this_well_zmws <- lima_zmws[are_this_well]
    fasta_folder <- file.path(fasta_output_dir, well_names[[i]])
    fastq_folder <- file.path(fastq_output_dir, well_names[[i]])
    dir.create(fasta_folder, showWarnings = FALSE)
    dir.create(fastq_folder, showWarnings = FALSE)
    if (length(this_well_zmws) > max_num_ZMWs) {
      this_well_zmws <- this_well_zmws[seq_len(max_num_ZMWs)]
    }
    for (zmw in this_well_zmws) {
      # message(paste0("Exporting reads for ZMW ", zmw))
      ExportReadsForZMW(zmw, output_dir = fasta_folder, use_fastq = FALSE)
      ExportReadsForZMW(zmw, output_dir = fastq_folder, use_fastq = TRUE)
    }
  }
  message("")
}



