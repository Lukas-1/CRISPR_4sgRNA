### 29th July 2021 ###



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

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(s2r2_R_objects_directory, "05) Read in PacBio data.RData"))
# load(file.path(s2r2_R_objects_directory, "06) Perform pairwise alignments with the reference sequence.RData"))




# Extract barcodes --------------------------------------------------------

are_eligible <- (ccs_df[["Well_exists"]] %in% TRUE) &
                (ccs_df[["Read_quality"]] > 0)

dummy_columns <-  c("ZMW", "Combined_ID", "Plate_number", "Well_number",
                    "Well_barcodes_orientation_fwd"
                    )
dummy_alignments_df <- ccs_df[are_eligible, dummy_columns]
row.names(dummy_alignments_df) <- NULL
names(dummy_alignments_df)[names(dummy_alignments_df) == "Well_barcodes_orientation_fwd"] <- "Orientation_fwd"

barcodes_df <- GetBarcodes(ccs_df, dummy_alignments_df)





# Save data ---------------------------------------------------------------

save(list = "barcodes_df",
     file = file.path(s2r2_R_objects_directory, "07) Extract barcode sequences and quality scores.RData")
     )




