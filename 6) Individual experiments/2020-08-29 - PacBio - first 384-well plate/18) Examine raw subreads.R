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

load(file.path(R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(R_objects_directory, "05) Read in PacBio data - consensus reads - ccs3.RData"))
load(file.path(R_objects_directory, "05) Read in PacBio data - demultiplexed - ccs3.RData"))
load(file.path(R_objects_directory, "05) Read in PacBio data - ccs5 ZMWs.RData"))
load(file.path(R_objects_directory, "17) Import raw subreads.RData"))





# Export subreads (before consensus calling) ------------------------------

ExportReadsForZMW(4522423)

message("Exporting reads for SmrtLink7_CCS3...")
ExportSubreadsForWells(fasta_output_dir = file.path(subreads_fasta_directory, "SmrtLink7_CCS3"),
                       fastq_output_dir = file.path(subreads_fastq_directory, "SmrtLink7_CCS3")
                       )












