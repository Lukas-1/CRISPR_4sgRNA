### 31st December 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "01) Define titles and labels.R"))
source(file.path(R_functions_directory, "16) Showing metrics in a schematic of a 384-well plate.R"))




# Define folder paths -----------------------------------------------------

s2rI_directory           <- file.path(experiments_directory, "2021-12-08 - integrate PacBio data")
p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")
s2rI_R_objects_directory <- file.path(s2rI_directory, "3) R objects")
file_output_directory    <- file.path(s2rI_directory, "5) Output")
plots_output_directory   <- file.path(file_output_directory, "Figures", "Schemas")
PNGs_output_directory    <- file.path(file_output_directory, "PNGs", "Schemas")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(s2rI_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2rI_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(s2rI_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))




# General preparation -----------------------------------------------------

sg_sequences_df[["Empty_well"]] <- FALSE




# Export all plates -------------------------------------------------------

use_plate_numbers <- plates_df[["Plate_number"]]#[order(plates_df[["Plate_rank"]])]

DrawSchematicsForAllPlates(export_PNGs = TRUE)






