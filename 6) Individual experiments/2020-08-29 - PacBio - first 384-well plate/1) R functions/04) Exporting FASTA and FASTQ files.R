### 20th October 2020 ###



# Define functions --------------------------------------------------------

ExportSequences <- function(lima_reads,
                            report_df,
                            fasta_output_dir,
                            fastq_output_dir,
                            ccs_reads = NULL,
                            append_to_file_name = "",
                            prefer_ccs = TRUE,
                            use_zmws = NULL
                            ) {

  stopifnot("sg_sequences_df" %in% ls(envir = globalenv()))

  file_names <- paste0("well",
                       formatC(seq_len(384), flag = "0", width = 3),
                       append_to_file_name
                       )

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

  for (i in sg_sequences_df[["Well_number"]]) {
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

    fasta_path <- file.path(fasta_output_dir, paste0(file_names[[i]], ".fasta"))
    fastq_path <- file.path(fastq_output_dir, paste0(file_names[[i]], ".fastq"))

    writeXStringSet(export_seq, filepath = fasta_path)
    writeQualityScaledXStringSet(export_fastq, filepath = fastq_path)
  }
}
