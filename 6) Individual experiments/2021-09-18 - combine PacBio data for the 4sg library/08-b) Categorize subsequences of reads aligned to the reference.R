### 27th September 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "07) Categorizing subsequences of reads aligned to the reference.R"))




# Define folder paths -----------------------------------------------------

s2rC_directory           <- file.path(experiments_directory, "2021-09-18 - combine PacBio data for the 4sg library")
s2rC_R_objects_directory <- file.path(s2rC_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(s2rC_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(s2rC_R_objects_directory, "08-a) Categorize subsequences of reads aligned to the reference.RData"))





# Extract sequences -------------------------------------------------------

extracted_df <- ProcessExtractedDf(extracted_df, unique_IDs = sg_sequences_df[["Combined_ID"]])




# Save data ---------------------------------------------------------------

save(list = c("extracted_df"),
     file = file.path(s2rC_R_objects_directory, "08-b) Categorize subsequences of reads aligned to the reference.RData")
     )




