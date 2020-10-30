### 27th October 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))
source(file.path(R_functions_directory, "12) Examining and exporting raw subreads.R"))




# Define folder paths -----------------------------------------------------

plate2_directory          <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-18 - PacBio - second 384-well plate")
p1_R_objects_directory    <- file.path(plate1_directory, "3) R objects")
p2_R_objects_directory    <- file.path(plate2_directory, "2) R objects")
file_output_directory     <- file.path(plate2_directory, "3) Output")
subreads_output_directory <- file.path(file_output_directory, "Subreads")
subreads_fasta_directory  <- file.path(subreads_output_directory, "Fasta")
subreads_fastq_directory  <- file.path(subreads_output_directory, "Fastq")
subreads_zmws_directory   <- file.path(subreads_output_directory, "Individual ZMWs")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(p2_R_objects_directory, "01) Import and process sgRNA sequences.RData"))
# load(file.path(p2_R_objects_directory, "03) Read in PacBio data - consensus reads - ccs3.RData"))
# load(file.path(p2_R_objects_directory, "03) Read in PacBio data - demultiplexed - ccs3.RData"))
# load(file.path(p2_R_objects_directory, "03) Read in PacBio data - ccs5 ZMWs.RData"))
load(file.path(p2_R_objects_directory, "03) Read in PacBio data - consensus reads - ccs5.RData"))
load(file.path(p2_R_objects_directory, "03) Read in PacBio data - demultiplexed - ccs5.RData"))
load(file.path(p2_R_objects_directory, "15) Import raw subreads.RData"))





# Export subreads (before consensus calling) ------------------------------

############################
### DELETE THIS LATER!!! ###
############################
sl7_ccs3_ccs <- sl7_ccs5_ccs
sl7_ccs3_lima <- sl7_ccs5_lima
sl7_ccs3_report_df <- sl7_ccs5_report_df
############################
############################
############################

ExportReadsForZMW(4588372)

message("Exporting reads for SmrtLink7_CCS3...")
ExportSubreadsForWells(fasta_output_dir = file.path(subreads_fasta_directory, "SmrtLink7_CCS5"),
                       fastq_output_dir = file.path(subreads_fastq_directory, "SmrtLink7_CCS5"),
                       wells_vec = sg_sequences_df[["Well_number"]]
                       )












