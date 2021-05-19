### 27th October 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
file_directory        <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))
source(file.path(R_functions_directory, "12) Examining and exporting raw subreads.R"))






# Define folder paths -----------------------------------------------------

R_objects_directory       <- file.path(file_directory, "3) R objects")
file_output_directory     <- file.path(file_directory, "5) Output")
subreads_output_directory <- file.path(file_output_directory, "Subreads")
subreads_fasta_directory  <- file.path(subreads_output_directory, "Fasta")
subreads_fastq_directory  <- file.path(subreads_output_directory, "Fastq")
subreads_zmws_directory   <- file.path(subreads_output_directory, "Individual ZMWs")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "05) Read in PacBio data.RData"))
load(file.path(R_objects_directory, "18) Import raw subreads.RData"))




# Export selected subreads ------------------------------------------------

ExportReadsForZMW(57672593, sl7_ccs_df)
ExportReadsForZMW(4522423, sl7_ccs_df)
ExportReadsForZMW(32965024, sl7_ccs_df)




# Export a small sample of raw subreads (before consensus calling) --------

message("Exporting reads for SmrtLink7_CCS3...")
ExportSubreadsForWells(sl7_ccs_df,
                       fasta_output_dir = file.path(subreads_fasta_directory, "SmrtLink7_CCS3"),
                       fastq_output_dir = file.path(subreads_fastq_directory, "SmrtLink7_CCS3"),
                       max_num_ZMWs = 20L
                       )












