## 2022-04-06


# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")

plate1_directory <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
p1_R_functions_directory <- file.path(plate1_directory, "1) R functions")
source(file.path(p1_R_functions_directory, "13) Importing PacBio reads.R"))



# Define folder paths -----------------------------------------------------

project_dir <- file.path(experiments_directory, "2022-04-06 - PacBio pooled 4sg - first trial")
rdata_dir   <- file.path(project_dir, "03_R_objects")
input_dir   <- file.path(project_dir, "02_input_data")
sam_path    <- file.path(input_dir, "m64141e_220403_204315.reads.sam")



# Read in data ------------------------------------------------------------

ccs_sam_file <- ProcessSAM(sam_path)




# Filter data -------------------------------------------------------------

are_not_minus_one <- ccs_sam_file[, "Read_quality"] != -1
ccs_df <- ccs_sam_file[are_not_minus_one, 1:5]
row.names(ccs_df) <- NULL




# Save data ---------------------------------------------------------------

save(ccs_df, file = file.path(rdata_dir, "01_read_in_data.RData"))



