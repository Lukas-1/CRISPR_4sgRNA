### 20th October 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))
source(file.path(R_functions_directory, "06) Extracting barcodes.R"))



# Define folder paths -----------------------------------------------------

plate2_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-18 - PacBio - second 384-well plate")
p1_R_objects_directory <- file.path(plate1_directory, "3) R objects")
p2_R_objects_directory <- file.path(plate2_directory, "2) R objects")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(p2_R_objects_directory, "02) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(p2_R_objects_directory, "03) Read in PacBio data.RData"))
load(file.path(p2_R_objects_directory, "04) Perform pairwise alignments with the reference sequence.RData"))




# Extract barcodes --------------------------------------------------------

sl7_barcodes_df <- GetBarcodes(sl7_ccs_df, sl7_alignments_df, wells_vec = sg_sequences_df[["Well_number"]])
sl9_barcodes_df <- GetBarcodes(sl9_ccs_df, sl9_alignments_df, wells_vec = sg_sequences_df[["Well_number"]])





# Save data ---------------------------------------------------------------

save(list = paste0("sl", c(7, 9), "_barcodes_df"),
     file = file.path(p2_R_objects_directory, "05) Extract barcode sequences and quality scores.RData")
     )












