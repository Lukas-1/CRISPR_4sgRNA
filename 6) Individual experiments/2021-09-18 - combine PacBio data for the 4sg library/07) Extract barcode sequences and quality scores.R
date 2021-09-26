### 24th September 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory      <- "~/CRISPR"
experiments_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory           <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
s2r1_directory             <- file.path(experiments_directory, "2021-04-03 - PacBio - first Sequel-II run")
p1_R_functions_directory   <- file.path(plate1_directory, "1) R functions")
s2r1_R_functions_directory <- file.path(s2r1_directory, "1) R functions")

source(file.path(p1_R_functions_directory, "02) Analyzing reads.R"))
source(file.path(s2r1_R_functions_directory, "02) Extracting barcodes.R"))




# Define folder paths -----------------------------------------------------

s2rC_directory           <- file.path(experiments_directory, "2021-09-18 - combine PacBio data for the 4sg library")
p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")
s2rC_R_objects_directory <- file.path(s2rC_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(s2rC_R_objects_directory, "05) Read in PacBio data.RData"))
load(file.path(s2rC_R_objects_directory, "06) Perform pairwise alignments with the reference sequence.RData"))




# Extract barcodes --------------------------------------------------------

barcodes_df <- GetBarcodes(ccs_df,
                           alignments_df,
                           tolerate_unexpected_barcodes = TRUE
                           )



# Save data ---------------------------------------------------------------

save(list = "barcodes_df",
     file = file.path(s2rC_R_objects_directory, "07) Extract barcode sequences and quality scores.RData")
     )




