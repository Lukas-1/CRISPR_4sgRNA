### 19th September 2020 ###




# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "01) Define titles and labels.R"))
source(file.path(R_functions_directory, "16) Showing metrics in a schematic of a 384-well plate.R"))




# Define folder paths -----------------------------------------------------

plate2_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-18 - PacBio - second 384-well plate")
p1_R_objects_directory <- file.path(plate1_directory, "3) R objects")
p2_R_objects_directory <- file.path(plate2_directory, "2) R objects")
file_output_directory  <- file.path(plate2_directory, "3) Output")
plots_output_directory <- file.path(file_output_directory, "Figures")





# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(p2_R_objects_directory, "01) Import and process sgRNA sequences.RData"))
load(file.path(p2_R_objects_directory, "09) Process demultiplexed PacBio reads.RData"))





# Exclude problematic wells, for the time being ---------------------------

sg_sequences_df[["Empty_well"]] <- ifelse(sg_sequences_df[["Well_number"]] == 2,
                                          TRUE, FALSE
                                          )




# Display metrics in the layout of a 384-well plate -----------------------

DrawAllSchematicsForOnePlate()

BarPlotPanel(sl7_ccs3_df_list[["original_summary_df"]],
             "Longest_subsequence",
             sg_sequences_df,
             indicate_homologies = TRUE
             )





