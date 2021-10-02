### 30th September 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R")) # For GetWellNumbers
source(file.path(R_functions_directory, "04) Exporting FASTA and FASTQ files.R"))
source(file.path(R_functions_directory, "08) Processing demultiplexed PacBio reads.R")) # For GetFeaturesData



# Define folder paths -----------------------------------------------------

s2rC_directory           <- file.path(experiments_directory, "2021-09-18 - combine PacBio data")
s2rC_R_objects_directory <- file.path(s2rC_directory, "3) R objects")
file_output_directory    <- file.path(s2rC_directory, "5) Output")
fastq_output_directory   <- file.path(file_output_directory, "Fastq")




# Load data ---------------------------------------------------------------

load(file.path(s2rC_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(s2rC_R_objects_directory, "05) Read in PacBio data.RData"))
load(file.path(s2rC_R_objects_directory, "08) Categorize subsequences of reads aligned to the reference.RData"))
load(file.path(s2rC_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))





# Prepare for exporting sequences -----------------------------------------

ccs_df[["Category_string"]] <- MakeCategoryString(ccs_df, extracted_df)

ccs_df[["ZMW_string"]] <- paste0("S2R3P3_",
                                 ccs_df[["ZMW"]], "_",
                                 ccs_df[["Category_string"]]
                                 )

ccs_df[["Passed_filters"]] <- ccs_df[["Plate_passed_filters"]] &
                              (ccs_df[["Well_passed_filters"]] %in% TRUE)

ccs3_zmws <- GetCCS3_ZMWs(ccs_df)
ccs5_zmws <- GetCCS5_ZMWs(ccs_df)
ccs7_zmws <- GetCCS7_ZMWs(ccs_df)

reads_df <- ccs3_df_list[["individual_reads_df"]]

pass_bc   <- reads_df[["Passes_barcode_filters"]] == 1
pass_read <- reads_df[["Passes_read_quality"]] == 1
pass_sg   <- reads_df[["Passes_sg_quality"]] == 1

passing_bc_zmws   <- reads_df[["ZMW"]][pass_bc]
passing_read_zmws <- reads_df[["ZMW"]][pass_bc & pass_read]
# passing_sg_zmws   <- reads_df[["ZMW"]][pass_bc & pass_read & pass_sg]



# Export sequences --------------------------------------------------------

for (filter_reads in c("Unfiltered", "Filtered reads", "Filtered cross-plate contaminations")) { #, "Filtered reads" # "Filtered gRNAs"
  for (split_reads in c(FALSE, TRUE)) {

    if (filter_reads == "Unfiltered") {
      first_half <- "a) Unfiltered"
      use_ccs3_zmws <- ccs3_zmws
      use_ccs5_zmws <- ccs5_zmws
      use_ccs7_zmws <- ccs7_zmws
    } else if (filter_reads == "Filtered reads") {
      first_half <- "b) Filtered reads"
      use_ccs3_zmws <- passing_read_zmws
      use_ccs5_zmws <- intersect(ccs5_zmws, passing_read_zmws)
      use_ccs7_zmws <- intersect(ccs7_zmws, passing_read_zmws)
    } else if (filter_reads == "Filtered cross-plate contaminations") {
      first_half <- "c) Filtered cross-plate"
      use_ccs3_zmws <- passing_read_zmws
      use_ccs5_zmws <- intersect(ccs5_zmws, passing_read_zmws)
      use_ccs7_zmws <- intersect(ccs7_zmws, passing_read_zmws)
    }

    if (split_reads) {
      sub_folder <- paste0(first_half, " - split into chunks")
    } else {
      sub_folder <- paste0(first_half, " - all reads")
    }

    message(paste0("Exporting reads into the '", sub_folder, "' folders..."))
    message("Exporting reads for CCS3...")
    ExportSequences(ccs_df,
                    fasta_output_dir    = NULL,
                    fastq_output_dir    = file.path(fastq_output_directory, sub_folder, "CCS3"),
                    append_to_file_name = "_ccs3",
                    use_zmws            = use_ccs3_zmws,
                    split_into_chunks   = split_reads,
                    ID_column           = "Combined_ID",
                    unique_IDs          = library_df[["Combined_ID"]],
                    export_fasta        = FALSE
                    )

    message("Exporting reads for CCS7...")
    ExportSequences(ccs_df,
                    fasta_output_dir    = NULL,
                    fastq_output_dir    = file.path(fastq_output_directory, sub_folder, "CCS7"),
                    append_to_file_name = "_ccs7",
                    use_zmws            = use_ccs7_zmws,
                    split_into_chunks   = split_reads,
                    ID_column           = "Combined_ID",
                    unique_IDs          = library_df[["Combined_ID"]],
                    export_fasta        = FALSE
                    )

  }
}



