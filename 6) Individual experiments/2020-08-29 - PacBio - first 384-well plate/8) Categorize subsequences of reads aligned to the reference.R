### 15th October 2020 ###



# Import packages and source code -----------------------------------------

library("Rsamtools")

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "2) Functions for analyzing reads.R"))




# Define folder paths -----------------------------------------------------

R_objects_directory <- file.path(file_directory, "3) R objects")



# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "1) Process and export barcodes.RData"))
load(file.path(R_objects_directory, "4) Create reference sequences for each well - raw sequences.RData"))
load(file.path(R_objects_directory, "5) Read in PacBio data - consensus reads.RData"))
load(file.path(R_objects_directory, "5) Read in PacBio data - demultiplexed.RData"))
load(file.path(R_objects_directory, "7) Process demultiplexed PacBio reads.RData"))





# Define functions --------------------------------------------------------




# Extract barcodes --------------------------------------------------------

sl7_ccs3_barcodes_df <- GetBarcodes(use_sl7 = TRUE, use_ccs3 = TRUE)
sl7_ccs5_barcodes_df <- GetBarcodes(use_sl7 = TRUE, use_ccs3 = FALSE)
sl9_ccs3_barcodes_df <- GetBarcodes(use_sl7 = FALSE, use_ccs3 = TRUE)
sl9_ccs5_barcodes_df <- GetBarcodes(use_sl7 = FALSE, use_ccs3 = FALSE)






# Save data ---------------------------------------------------------------

save(list = paste0(c("sl7_ccs3", "sl7_ccs5", "sl9_ccs3", "sl9_ccs5"), "_barcodes_df"),
     file = file.path(R_objects_directory, "7) Extract barcode sequences and quality scores.RData")
     )











