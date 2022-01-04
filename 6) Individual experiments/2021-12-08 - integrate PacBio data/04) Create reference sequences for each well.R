### 8th December 2021 ###




# Import packages and source code -----------------------------------------




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")

s2rC_directory           <- file.path(experiments_directory, "2021-09-18 - combine PacBio data")
s2rI_directory           <- file.path(experiments_directory, "2021-12-08 - integrate PacBio data")

s2rC_R_objects_directory <- file.path(s2rC_directory, "3) R objects")
s2rI_R_objects_directory <- file.path(s2rI_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(s2rC_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))




# Save data ---------------------------------------------------------------

save(list = "sg_sequences_df",
     file = file.path(s2rI_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData")
     )



