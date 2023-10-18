### 15th August 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory       <- "~/CRISPR_4sgRNA"
general_functions_directory <- file.path(CRISPR_root_directory, "1) R scripts/1) R functions")
experiments_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory            <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory       <- file.path(plate1_directory, "1) R functions")

source(file.path(general_functions_directory, "05) Mapping sequences to the human genome.R"))
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))
source(file.path(general_functions_directory, "31) Finding 0MM hits for variable-length sgRNAs.R"))

source(file.path(R_functions_directory, "24) Finding unintended targets of mutated gRNAs.R"))




# Define folder paths -----------------------------------------------------

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(s2r2_R_objects_directory, "08) Categorize subsequences of reads aligned to the reference.RData"))
load(file.path(s2r2_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))




# Define all "mutated" sequences to be considered -------------------------

mutations_df <- FilterMutations(extracted_df, ccs7_df_list)





# Annotate mutated sequences with 0MM off-target sites --------------------

mutations_df <- Find19or20MutatedAndTemplate(mutations_df)





# Save data ---------------------------------------------------------------

save(list = "mutations_df",
     file = file.path(s2r2_R_objects_directory, "26) Annotate mutated sgRNAs with any perfect matches in the genome.RData")
     )



