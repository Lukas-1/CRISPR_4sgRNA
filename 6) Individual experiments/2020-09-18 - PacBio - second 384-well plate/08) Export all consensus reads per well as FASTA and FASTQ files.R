### 20th October 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory  <- "~/CRISPR"
plate1_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R")) # For GetWellNumbers
source(file.path(R_functions_directory, "04) Exporting FASTA and FASTQ files.R"))






# Define folder paths -----------------------------------------------------

plate2_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-18 - PacBio - second 384-well plate")
p1_R_objects_directory <- file.path(plate1_directory, "3) R objects")
p2_R_objects_directory <- file.path(plate2_directory, "2) R objects")

file_output_directory  <- file.path(plate2_directory, "3) Output")
fasta_output_directory <- file.path(file_output_directory, "Fasta")
fastq_output_directory <- file.path(file_output_directory, "Fastq")





# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(p2_R_objects_directory, "02) Create reference sequences for each well - raw sequences.RData"))
# load(file.path(p2_R_objects_directory, "03) Read in PacBio data - demultiplexed - ccs3.RData"))
# load(file.path(p2_R_objects_directory, "03) Read in PacBio data - consensus reads - ccs3.RData"))
load(file.path(p2_R_objects_directory, "03) Read in PacBio data - demultiplexed - ccs5.RData"))
load(file.path(p2_R_objects_directory, "03) Read in PacBio data - consensus reads - ccs5.RData"))
# load(file.path(p2_R_objects_directory, "03) Read in PacBio data - ccs5 ZMWs.RData"))
load(file.path(p2_R_objects_directory, "07) Process demultiplexed PacBio reads.RData"))





# Prepare for exporting sequences -----------------------------------------

# sl7_reads_df <- sl7_ccs3_df_list[["individual_reads_df"]]
# sl9_reads_df <- sl9_ccs3_df_list[["individual_reads_df"]]
sl7_reads_df <- sl7_ccs5_df_list[["individual_reads_df"]]

sl7_passing_zmws <- sl7_reads_df[["ZMW"]][sl7_reads_df[["Passes_barcode_filters"]] == 1]
# sl9_passing_zmws <- sl9_reads_df[["ZMW"]][sl9_reads_df[["Passes_barcode_filters"]] == 1]




# Export sequences --------------------------------------------------------

for (filter_reads in c(FALSE, TRUE)) {
  for (split_reads in c(FALSE, TRUE)) {

    if (filter_reads) {
      first_half <- "Filtered"
      # use_sl7_ccs3_zmws <- sl7_passing_zmws
      # use_sl7_ccs5_zmws <- intersect(sl7_ccs5_lima_zmws, sl7_passing_zmws)
      # use_sl9_ccs3_zmws <- sl9_passing_zmws
      # use_sl9_ccs5_zmws <- intersect(sl9_ccs5_lima_zmws, sl9_passing_zmws)
      use_sl7_ccs5_zmws <- sl7_passing_zmws
    } else {
      first_half <- "Unfiltered"
      # use_sl7_ccs3_zmws <- NULL
      # use_sl7_ccs5_zmws <- sl7_ccs5_lima_zmws
      # use_sl9_ccs3_zmws <- NULL
      # use_sl9_ccs5_zmws <- sl9_ccs5_lima_zmws
      use_sl7_ccs5_zmws <- NULL
    }

    if (split_reads) {
      sub_folder <- paste0(first_half, "_split_into_chunks")
    } else {
      sub_folder <- paste0(first_half, "_all_reads")
    }

    message(paste0("Exporting reads into the '", sub_folder, "' subfolders..."))

    ExportSequences(sl7_ccs5_lima,
                    sl7_ccs5_report_df,
                    fasta_output_dir = file.path(fasta_output_directory, sub_folder, "SmrtLink7_CCS5"),
                    fastq_output_dir = file.path(fastq_output_directory, sub_folder, "SmrtLink7_CCS5"),
                    append_to_file_name = "_ccs5",
                    ccs_reads = sl7_ccs5_ccs,
                    use_zmws = use_sl7_ccs5_zmws,
                    split_into_chunks = split_reads,
                    wells_vec = sg_sequences_df[["Well_number"]]
                    )


    # message("Exporting reads for SmrtLink7_CCS3...")
    # ExportSequences(sl7_ccs3_lima,
    #                 sl7_ccs3_report_df,
    #                 fasta_output_dir = file.path(fasta_output_directory, sub_folder, "SmrtLink7_CCS3"),
    #                 fastq_output_dir = file.path(fastq_output_directory, sub_folder, "SmrtLink7_CCS3"),
    #                 append_to_file_name = "_ccs3",
    #                 ccs_reads = sl7_ccs3_ccs,
    #                 use_zmws = use_sl7_ccs3_zmws,
    #                 split_into_chunks = split_reads,
    #                 wells_vec = sg_sequences_df[["Well_number"]]
    #                 )
    #
    # message("Exporting reads for SmrtLink7_CCS5...")
    # ExportSequences(sl7_ccs3_lima,
    #                 sl7_ccs3_report_df,
    #                 fasta_output_dir = file.path(fasta_output_directory, sub_folder, "SmrtLink7_CCS5"),
    #                 fastq_output_dir = file.path(fastq_output_directory, sub_folder, "SmrtLink7_CCS5"),
    #                 append_to_file_name = "_ccs5",
    #                 ccs_reads = sl7_ccs3_ccs,
    #                 use_zmws = use_sl7_ccs5_zmws,
    #                 split_into_chunks = split_reads,
    #                 wells_vec = sg_sequences_df[["Well_number"]]
    #                 )
    #
    # message("Exporting reads for SmrtLink9_CCS3...")
    # ExportSequences(sl9_ccs3_lima,
    #                 sl9_ccs3_report_df,
    #                 fasta_output_dir = file.path(fasta_output_directory, sub_folder, "SmrtLink9_CCS3"),
    #                 fastq_output_dir = file.path(fastq_output_directory, sub_folder, "SmrtLink9_CCS3"),
    #                 append_to_file_name = "_ccs3",
    #                 ccs_reads = sl9_ccs3_ccs,
    #                 use_zmws = use_sl9_ccs3_zmws,
    #                 split_into_chunks = split_reads,
    #                 wells_vec = sg_sequences_df[["Well_number"]]
    #                 )
    #
    # message("Exporting reads for SmrtLink9_CCS5...")
    # ExportSequences(sl9_ccs3_lima,
    #                 sl9_ccs3_report_df,
    #                 fasta_output_dir = file.path(fasta_output_directory, sub_folder, "SmrtLink9_CCS5"),
    #                 fastq_output_dir = file.path(fastq_output_directory, sub_folder, "SmrtLink9_CCS5"),
    #                 append_to_file_name = "_ccs5",
    #                 ccs_reads = sl9_ccs3_ccs,
    #                 use_zmws = use_sl9_ccs5_zmws,
    #                 split_into_chunks = split_reads,
    #                 wells_vec = sg_sequences_df[["Well_number"]]
    #                 )

  }
}



