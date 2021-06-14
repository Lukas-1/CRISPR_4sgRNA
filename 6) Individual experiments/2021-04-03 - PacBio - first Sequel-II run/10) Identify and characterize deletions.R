### 1st June 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "07) Categorizing subsequences of reads aligned to the reference.R"))
source(file.path(R_functions_directory, "17) Characterizing deletions.R"))




# Define folder paths -----------------------------------------------------

sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(sql2_R_objects_directory, "06) Perform pairwise alignments with the reference sequence.RData"))





# Identify and characterize large deletions -------------------------------

deletions_df <- CompileDeletions(alignments_df,
                                 "Combined_ID",
                                 sg_sequences_df[["Combined_ID"]]
                                 )




# Save data ---------------------------------------------------------------

save(list = "deletions_df",
     file = file.path(sql2_R_objects_directory, "10) Identify and characterize deletions.RData")
     )




