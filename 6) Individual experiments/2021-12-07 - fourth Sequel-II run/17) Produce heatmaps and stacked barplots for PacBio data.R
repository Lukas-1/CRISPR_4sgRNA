### 13th December 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "09) Producing heatmaps.R"))
source(file.path(R_functions_directory, "11) Creating stacked barplots for visualizing alterations.R"))




# Define folder paths -----------------------------------------------------

s2r4_directory           <- file.path(experiments_directory, "2021-12-07 - fourth Sequel-II run")
s2r4_R_objects_directory <- file.path(s2r4_directory, "3) R objects")
file_output_directory    <- file.path(s2r4_directory, "5) Output")
plots_output_directory   <- file.path(file_output_directory, "Figures")
# PNGs_output_directory    <- file.path(file_output_directory, "PNGs")




# Load data ---------------------------------------------------------------

load(file.path(s2r4_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2r4_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(s2r4_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))





# Export individual graphics ----------------------------------------------

use_plate_numbers <- plates_df[["Plate_number"]]

DrawBarplotsAndHeatmapsForAllPlates(export_PNGs = FALSE)







