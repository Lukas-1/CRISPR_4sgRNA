### 13th October 2020 ###




# Import packages and source code -----------------------------------------

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R")) # For GetCCS5_ZMWs
source(file.path(R_functions_directory, "04) Exporting FASTA and FASTQ files.R"))
source(file.path(R_functions_directory, "08) Processing demultiplexed PacBio reads.R")) # For GetFeaturesData





# Define folder paths -----------------------------------------------------

R_objects_directory    <- file.path(file_directory, "3) R objects")
file_output_directory  <- file.path(file_directory, "5) Output")
fasta_output_directory <- file.path(file_output_directory, "Fasta")
fastq_output_directory <- file.path(file_output_directory, "Fastq")





# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(R_objects_directory, "05) Read in PacBio data.RData"))
load(file.path(R_objects_directory, "08) Categorize subsequences of reads aligned to the reference.RData"))
load(file.path(R_objects_directory, "09) Process demultiplexed PacBio reads.RData"))





# Prepare for exporting sequences -----------------------------------------

sl7_ccs_df[["Category_string"]] <- MakeCategoryString(sl7_ccs_df, sl7_extracted_df)
sl9_ccs_df[["Category_string"]] <- MakeCategoryString(sl9_ccs_df, sl9_extracted_df)

sl7_ccs5_lima_zmws <- GetCCS5_ZMWs(sl7_ccs_df)
sl9_ccs5_lima_zmws <- GetCCS5_ZMWs(sl9_ccs_df)

sl7_reads_df <- sl7_ccs3_df_list[["individual_reads_df"]]
sl9_reads_df <- sl9_ccs3_df_list[["individual_reads_df"]]

sl7_pass_bc <- sl7_reads_df[["Passes_barcode_filters"]] == 1
sl9_pass_bc <- sl9_reads_df[["Passes_barcode_filters"]] == 1

sl7_pass_read <- sl7_reads_df[["Passes_read_quality"]] == 1
sl9_pass_read <- sl9_reads_df[["Passes_read_quality"]] == 1

sl7_pass_sg <- sl7_reads_df[["Passes_sg_quality"]] == 1
sl9_pass_sg <- sl9_reads_df[["Passes_sg_quality"]] == 1

sl7_passing_bc_zmws <- sl7_reads_df[["ZMW"]][sl7_pass_bc]
sl9_passing_bc_zmws <- sl9_reads_df[["ZMW"]][sl9_pass_bc]

sl7_passing_read_zmws <- sl7_reads_df[["ZMW"]][sl7_pass_bc & sl7_pass_read]
sl9_passing_read_zmws <- sl9_reads_df[["ZMW"]][sl9_pass_bc & sl9_pass_read]

sl7_passing_sg_zmws <- sl7_reads_df[["ZMW"]][sl7_pass_bc & sl7_pass_read & sl7_pass_sg]
sl9_passing_sg_zmws <- sl9_reads_df[["ZMW"]][sl9_pass_bc & sl9_pass_read & sl9_pass_sg]




# Export sequences --------------------------------------------------------

for (filter_reads in c("Unfiltered", "Filtered barcodes", "Filtered reads", "Filtered gRNAs")) {
  for (split_reads in c(FALSE, TRUE)) {

    if (filter_reads == "Unfiltered") {
      first_half <- "a) Unfiltered"
      use_sl7_ccs3_zmws <- NULL
      use_sl7_ccs5_zmws <- sl7_ccs5_lima_zmws
      use_sl9_ccs3_zmws <- NULL
      use_sl9_ccs5_zmws <- sl9_ccs5_lima_zmws
    } else if (filter_reads == "Filtered barcodes") {
      first_half <- "b) Filtered barcodes"
      use_sl7_ccs3_zmws <- sl7_passing_bc_zmws
      use_sl7_ccs5_zmws <- intersect(sl7_ccs5_lima_zmws, sl7_passing_bc_zmws)
      use_sl9_ccs3_zmws <- sl9_passing_bc_zmws
      use_sl9_ccs5_zmws <- intersect(sl9_ccs5_lima_zmws, sl9_passing_bc_zmws)
    } else if (filter_reads == "Filtered reads") {
      first_half <- "c) Filtered reads"
      use_sl7_ccs3_zmws <- sl7_passing_read_zmws
      use_sl7_ccs5_zmws <- intersect(sl7_ccs5_lima_zmws, sl7_passing_read_zmws)
      use_sl9_ccs3_zmws <- sl9_passing_read_zmws
      use_sl9_ccs5_zmws <- intersect(sl9_ccs5_lima_zmws, sl9_passing_read_zmws)
    } else if (filter_reads == "Filtered gRNAs") {
      first_half <- "d) Filtered gRNAs"
      use_sl7_ccs3_zmws <- sl7_passing_sg_zmws
      use_sl7_ccs5_zmws <- intersect(sl7_ccs5_lima_zmws, sl7_passing_sg_zmws)
      use_sl9_ccs3_zmws <- sl9_passing_sg_zmws
      use_sl9_ccs5_zmws <- intersect(sl9_ccs5_lima_zmws, sl9_passing_sg_zmws)
    }

    if (split_reads) {
      sub_folder <- paste0(first_half, " - split into chunks")
    } else {
      sub_folder <- paste0(first_half, " - all reads")
    }

    message(paste0("Exporting reads into the '", sub_folder, "' folders..."))


    message("Exporting reads for SmrtLink7_CCS3...")
    ExportSequences(sl7_ccs_df,
                    fasta_output_dir = file.path(fasta_output_directory, sub_folder, "SmrtLink7_CCS3"),
                    fastq_output_dir = file.path(fastq_output_directory, sub_folder, "SmrtLink7_CCS3"),
                    append_to_file_name = "_ccs3",
                    use_zmws = use_sl7_ccs3_zmws,
                    split_into_chunks = split_reads
                    )

    message("Exporting reads for SmrtLink7_CCS5...")
    ExportSequences(sl7_ccs_df,
                    fasta_output_dir = file.path(fasta_output_directory, sub_folder, "SmrtLink7_CCS5"),
                    fastq_output_dir = file.path(fastq_output_directory, sub_folder, "SmrtLink7_CCS5"),
                    append_to_file_name = "_ccs5",
                    use_zmws = use_sl7_ccs5_zmws,
                    split_into_chunks = split_reads
                    )

    message("Exporting reads for SmrtLink9_CCS3...")
    ExportSequences(sl9_ccs_df,
                    fasta_output_dir = file.path(fasta_output_directory, sub_folder, "SmrtLink9_CCS3"),
                    fastq_output_dir = file.path(fastq_output_directory, sub_folder, "SmrtLink9_CCS3"),
                    append_to_file_name = "_ccs3",
                    use_zmws = use_sl9_ccs3_zmws,
                    split_into_chunks = split_reads
                    )

    message("Exporting reads for SmrtLink9_CCS5...")
    ExportSequences(sl9_ccs_df,
                    fasta_output_dir = file.path(fasta_output_directory, sub_folder, "SmrtLink9_CCS5"),
                    fastq_output_dir = file.path(fastq_output_directory, sub_folder, "SmrtLink9_CCS5"),
                    append_to_file_name = "_ccs5",
                    use_zmws = use_sl9_ccs5_zmws,
                    split_into_chunks = split_reads
                    )

  }
}



