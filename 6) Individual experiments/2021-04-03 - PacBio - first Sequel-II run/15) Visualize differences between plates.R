### 28 May 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory      <- "~/CRISPR_4sgRNA"
plate1_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
sql2_directory             <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
p1_R_functions_directory   <- file.path(plate1_directory, "1) R functions")
sql2_R_functions_directory <- file.path(sql2_directory, "1) R functions")

source(file.path(p1_R_functions_directory, "01) Define titles and labels.R"))
source(file.path(sql2_R_functions_directory, "04) Visualizing differences between plates.R"))



# Define folder paths -----------------------------------------------------

sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")
file_output_directory    <- file.path(sql2_directory, "5) Output")
plots_output_directory   <- file.path(file_output_directory, "Figures", "Compare plates")
PNGs_output_directory    <- file.path(file_output_directory, "PNGs", "Compare plates")




# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(sql2_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))




# Prepare titles ----------------------------------------------------------

titles_list <- c(list("Count_total" = "Number of reads per well"),
                 titles_list
                 )




# Export graphics ---------------------------------------------------------

DrawAllPlateComparisons()




