### 15th October 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))
source(file.path(R_functions_directory, "05) Performing pairwise alignments with the reference sequence.R"))




# Define folder paths -----------------------------------------------------

R_objects_directory <- file.path(file_directory, "3) R objects")



# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(R_objects_directory, "04) Create reference sequences for each well - raw sequences.RData"))
load(file.path(R_objects_directory, "05) Read in PacBio data - consensus reads - ccs3.RData"))
load(file.path(R_objects_directory, "05) Read in PacBio data - demultiplexed - ccs3.RData")) # for the report_df




# Extract barcodes --------------------------------------------------------

sl7_alignments_df <- ExtractAlignedSequences(use_sl7 = TRUE)
sl9_alignments_df <- ExtractAlignedSequences(use_sl7 = FALSE)




# Save data ---------------------------------------------------------------

save(list = paste0("sl", c(7, 9), "_alignments_df"),
     file = file.path(R_objects_directory, "06) Perform pairwise alignments with the reference sequence.RData")
     )








