### 20th December 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "17) Characterizing contaminations (using aligned reads).R"))




# Define folder paths -----------------------------------------------------

s2rI_directory           <- file.path(experiments_directory, "2021-12-08 - integrate PacBio data")
s2rI_R_objects_directory <- file.path(s2rI_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(s2rI_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(s2rI_R_objects_directory, "08-b) Categorize subsequences of reads aligned to the reference.RData"))




# Characterize contaminations (using aligned reads) -----------------------

contam_df <- CharacterizeContaminations(extracted_df, sg_sequences_df)




# Save data ---------------------------------------------------------------

save(list = "contam_df",
     file = file.path(s2rI_R_objects_directory, "09) Characterize contaminations (using aligned reads).RData")
     )




