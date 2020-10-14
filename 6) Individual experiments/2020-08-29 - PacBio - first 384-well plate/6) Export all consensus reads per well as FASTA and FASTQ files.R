### 13th October 2020 ###





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
file_directory            <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_objects_directory       <- file.path(file_directory, "3) R objects")

file_output_directory     <- file.path(file_directory, "5) Output")
fasta_output_directory    <- file.path(file_output_directory, "Fasta")
fastq_output_directory    <- file.path(file_output_directory, "Fastq")





# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "1) Process and export barcodes.RData"))
load(file.path(R_objects_directory, "3) Import and process sgRNA sequences.RData"))
load(file.path(R_objects_directory, "4) Create reference sequences for each well - raw sequences.RData"))
load(file.path(R_objects_directory, "5) Read in PacBio data - demultiplexed.RData"))
load(file.path(R_objects_directory, "5) Read in PacBio data - consensus reads.RData"))





# Define functions --------------------------------------------------------

# This is function is currently also replicated in the "Extract barcode sequences..." script
GetWellNumbers <- function(lima_report_df) {
  combo_IDs_vec <- paste0(lima_report_df[["IdxLowestNamed"]], "--",
                          lima_report_df[["IdxHighestNamed"]]
                          )
  barcodes_to_wells_map[combo_IDs_vec]
}


ExportSequences <- function(lima_reads,
                            report_df,
                            fasta_output_dir,
                            fastq_output_dir,
                            ccs_reads = NULL,
                            append_to_file_name = "",
                            prefer_ccs = TRUE
                            ) {

  file_names <- paste0("well",
                       formatC(seq_len(384), flag = "0", width = 3),
                       append_to_file_name
                       )

  lima_zmws <- as.integer(substr(lima_reads[["qname"]], 22, nchar(lima_reads[["qname"]]) - 4))

  ccs_zmws <- as.integer(substr(report_df[["ZMW"]], 22, nchar(report_df[["ZMW"]])))

  ccs_well_numbers <- GetWellNumbers(report_df)
  lima_well_numbers <- ccs_well_numbers[match(lima_zmws, ccs_zmws)]

  use_ccs <- prefer_ccs && (!(is.null(ccs_reads)))
  if (!(is.null(ccs_reads))) {
    stopifnot(identical(ccs_zmws, as.integer(substr(ccs_reads[["qname"]], 22, nchar(ccs_reads[["qname"]]) - 4))))
  }

  for (i in seq_len(384)) {
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







# Export sequences --------------------------------------------------------

ExportSequences(sl7_ccs3_lima,
                sl7_ccs3_report_df,
                fasta_output_dir = file.path(fasta_output_directory, "SmrtLink7_CCS3"),
                fastq_output_dir = file.path(fastq_output_directory, "SmrtLink7_CCS3"),
                append_to_file_name = "_ccs3",
                ccs_reads = sl7_ccs3_ccs
                )

ExportSequences(sl7_ccs5_lima,
                sl7_ccs5_report_df,
                fasta_output_dir = file.path(fasta_output_directory, "SmrtLink7_CCS5"),
                fastq_output_dir = file.path(fastq_output_directory, "SmrtLink7_CCS5"),
                append_to_file_name = "_ccs5",
                ccs_reads = sl7_ccs5_ccs
                )

ExportSequences(sl9_ccs3_lima,
                sl9_ccs3_report_df,
                fasta_output_dir = file.path(fasta_output_directory, "SmrtLink9_CCS3"),
                fastq_output_dir = file.path(fastq_output_directory, "SmrtLink9_CCS3"),
                append_to_file_name = "_ccs3",
                ccs_reads = sl9_ccs3_ccs
                )

ExportSequences(sl9_ccs5_lima,
                sl9_ccs5_report_df,
                fasta_output_dir = file.path(fasta_output_directory, "SmrtLink9_CCS5"),
                fastq_output_dir = file.path(fastq_output_directory, "SmrtLink9_CCS5"),
                append_to_file_name = "_ccs5",
                ccs_reads = sl9_ccs5_ccs
                )
















