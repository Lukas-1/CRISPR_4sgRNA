### 14th June 2021 ###




# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "17) Characterizing contaminations (using aligned reads).R"))




# Define folder paths -----------------------------------------------------

plate2_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-18 - PacBio - second 384-well plate")
p2_R_objects_directory <- file.path(plate2_directory, "2) R objects")




# Load data ---------------------------------------------------------------

load(file.path(p2_R_objects_directory, "02) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(p2_R_objects_directory, "06) Categorize subsequences of reads aligned to the reference.RData"))





# Characterize contaminations (using aligned reads) -----------------------

sl7_contam_df <- CharacterizeContaminations(sl7_extracted_df, sg_sequences_df)
sl9_contam_df <- CharacterizeContaminations(sl9_extracted_df, sg_sequences_df)




# Save data ---------------------------------------------------------------

save(list = c("sl7_contam_df", "sl9_contam_df"),
     file = file.path(p2_R_objects_directory, "07) Characterize contaminations (using aligned reads).RData")
     )





