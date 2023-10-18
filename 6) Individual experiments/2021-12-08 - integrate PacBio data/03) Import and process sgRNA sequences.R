### 8th December 2021 ###




# Import packages and source code -----------------------------------------




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")

s2rC_directory           <- file.path(experiments_directory, "2021-09-18 - combine PacBio data")
s2rI_directory           <- file.path(experiments_directory, "2021-12-08 - integrate PacBio data")

s2rC_R_objects_directory <- file.path(s2rC_directory, "3) R objects")
s2rI_R_objects_directory <- file.path(s2rI_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(s2rC_R_objects_directory, "03) Import and process sgRNA sequences.RData"))





# Save data ---------------------------------------------------------------

save(list = "library_df",
     file = file.path(s2rI_R_objects_directory, "03) Import and process sgRNA sequences.RData")
     )



