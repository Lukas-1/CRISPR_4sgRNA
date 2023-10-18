### 1st June 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "07) Categorizing subsequences of reads aligned to the reference.R"))
source(file.path(R_functions_directory, "18) Characterizing deletions.R"))




# Define folder paths -----------------------------------------------------

plate2_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-18 - PacBio - second 384-well plate")
p2_R_objects_directory <- file.path(plate2_directory, "2) R objects")




# Load data ---------------------------------------------------------------

load(file.path(p2_R_objects_directory, "02) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(p2_R_objects_directory, "04) Perform pairwise alignments with the reference sequence.RData"))





# Identify and characterize deletions -------------------------------------

sl7_deletions_df <- CompileDeletions(sl7_alignments_df)
sl9_deletions_df <- CompileDeletions(sl9_alignments_df)




# Save data ---------------------------------------------------------------

save(list = c("sl7_deletions_df", "sl9_deletions_df"),
     file = file.path(p2_R_objects_directory, "08) Identify and characterize deletions.RData")
     )



