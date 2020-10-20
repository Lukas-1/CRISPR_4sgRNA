### 13th October 2020 ###




# Import packages and source code -----------------------------------------

library("Biostrings")

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R")) # For GetWellNumbers
source(file.path(R_functions_directory, "04) Exporting FASTA and FASTQ files.R"))






# Define folder paths -----------------------------------------------------

R_objects_directory    <- file.path(file_directory, "3) R objects")
file_output_directory  <- file.path(file_directory, "5) Output")
fasta_output_directory <- file.path(file_output_directory, "Fasta")
fastq_output_directory <- file.path(file_output_directory, "Fastq")





# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(R_objects_directory, "04) Create reference sequences for each well - raw sequences.RData"))
load(file.path(R_objects_directory, "05) Read in PacBio data - demultiplexed - ccs3.RData"))
load(file.path(R_objects_directory, "05) Read in PacBio data - consensus reads - ccs3.RData"))
load(file.path(R_objects_directory, "05) Read in PacBio data - ccs5 ZMWs.RData"))






# Export sequences --------------------------------------------------------

ExportSequences(sl7_ccs3_lima,
                sl7_ccs3_report_df,
                fasta_output_dir = file.path(fasta_output_directory, "SmrtLink7_CCS3"),
                fastq_output_dir = file.path(fastq_output_directory, "SmrtLink7_CCS3"),
                append_to_file_name = "_ccs3",
                ccs_reads = sl7_ccs3_ccs
                )

ExportSequences(sl7_ccs3_lima,
                sl7_ccs3_report_df,
                fasta_output_dir = file.path(fasta_output_directory, "SmrtLink7_CCS5"),
                fastq_output_dir = file.path(fastq_output_directory, "SmrtLink7_CCS5"),
                append_to_file_name = "_ccs5",
                ccs_reads = sl7_ccs3_ccs,
                use_zmws = sl7_ccs5_lima_zmws
                )

ExportSequences(sl9_ccs3_lima,
                sl9_ccs3_report_df,
                fasta_output_dir = file.path(fasta_output_directory, "SmrtLink9_CCS3"),
                fastq_output_dir = file.path(fastq_output_directory, "SmrtLink9_CCS3"),
                append_to_file_name = "_ccs3",
                ccs_reads = sl9_ccs3_ccs
                )

ExportSequences(sl9_ccs3_lima,
                sl9_ccs3_report_df,
                fasta_output_dir = file.path(fasta_output_directory, "SmrtLink9_CCS5"),
                fastq_output_dir = file.path(fastq_output_directory, "SmrtLink9_CCS5"),
                append_to_file_name = "_ccs5",
                ccs_reads = sl9_ccs3_ccs,
                use_zmws = sl9_ccs5_lima_zmws
                )











