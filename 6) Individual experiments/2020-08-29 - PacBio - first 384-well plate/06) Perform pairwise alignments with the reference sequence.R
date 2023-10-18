### 15th October 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
file_directory        <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "05) Performing pairwise alignments with the reference sequence.R"))




# Define folder paths -----------------------------------------------------

R_objects_directory <- file.path(file_directory, "3) R objects")



# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(R_objects_directory, "05) Read in PacBio data.RData"))




# Extract barcodes --------------------------------------------------------

sl7_alignments_df <- ExtractAlignedSequences(sl7_ccs_df, sg_sequences_df, parallel_mode = TRUE)
sl9_alignments_df <- ExtractAlignedSequences(sl9_ccs_df, sg_sequences_df, parallel_mode = TRUE)





# Save data ---------------------------------------------------------------

save(list = paste0("sl", c(7, 9), "_alignments_df"),
     file = file.path(R_objects_directory, "06) Perform pairwise alignments with the reference sequence.RData")
     )








