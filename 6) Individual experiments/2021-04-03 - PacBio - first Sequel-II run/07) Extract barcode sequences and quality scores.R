### 18th April 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory      <- "~/CRISPR"
plate1_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
p1_R_functions_directory   <- file.path(plate1_directory, "1) R functions")
sql2_directory             <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_functions_directory <- file.path(sql2_directory, "1) R functions")

source(file.path(p1_R_functions_directory, "02) Analyzing reads.R"))
source(file.path(sql2_R_functions_directory, "02) Extracting barcodes.R"))




# Define folder paths -----------------------------------------------------

p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(sql2_R_objects_directory, "05) Read in PacBio data.RData"))
load(file.path(sql2_R_objects_directory, "06) Perform pairwise alignments with the reference sequence.RData"))




# Extract barcodes --------------------------------------------------------

barcodes_df <- GetBarcodes(ccs_df, alignments_df)





# Save data ---------------------------------------------------------------

save(list = "barcodes_df",
     file = file.path(sql2_R_objects_directory, "07) Extract barcode sequences and quality scores.RData")
     )




