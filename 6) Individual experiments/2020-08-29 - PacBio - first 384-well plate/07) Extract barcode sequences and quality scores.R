### 15th October 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory  <- "~/CRISPR_4sgRNA"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))
source(file.path(R_functions_directory, "06) Extracting barcodes.R"))




# Define folder paths -----------------------------------------------------

R_objects_directory <- file.path(file_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(R_objects_directory, "05) Read in PacBio data.RData"))
load(file.path(R_objects_directory, "06) Perform pairwise alignments with the reference sequence.RData"))




# Extract barcodes --------------------------------------------------------

sl7_barcodes_df <- GetBarcodes(sl7_ccs_df, sl7_alignments_df)
sl9_barcodes_df <- GetBarcodes(sl9_ccs_df, sl9_alignments_df)





# Save data ---------------------------------------------------------------

save(list = paste0(c("sl7", "sl9"), "_barcodes_df"),
     file = file.path(R_objects_directory, "07) Extract barcode sequences and quality scores.RData")
     )




