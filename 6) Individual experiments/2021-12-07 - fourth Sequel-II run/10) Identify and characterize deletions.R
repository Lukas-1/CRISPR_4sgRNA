### 12th December 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")
source(file.path(R_functions_directory, "07) Categorizing subsequences of reads aligned to the reference.R"))
source(file.path(R_functions_directory, "18) Characterizing deletions.R"))




# Define folder paths -----------------------------------------------------

s2r4_directory           <- file.path(experiments_directory, "2021-12-07 - fourth Sequel-II run")
s2r4_R_objects_directory <- file.path(s2r4_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(s2r4_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(s2r4_R_objects_directory, "06) Perform pairwise alignments with the reference sequence.RData"))




# Identify and characterize large deletions -------------------------------

deletions_df <- CompileDeletions(alignments_df,
                                 "Combined_ID",
                                 sg_sequences_df[["Combined_ID"]]
                                 )



# Save data ---------------------------------------------------------------

save(list = "deletions_df",
     file = file.path(s2r4_R_objects_directory, "10) Identify and characterize deletions.RData")
     )




