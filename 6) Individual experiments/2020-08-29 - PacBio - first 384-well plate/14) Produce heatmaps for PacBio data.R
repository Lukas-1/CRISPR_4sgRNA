### 28 May 2020 ###




# Import packages and source code -----------------------------------------

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "09) Producing heatmaps.R"))




# Define folder paths -----------------------------------------------------

R_objects_directory    <- file.path(file_directory, "3) R objects")
file_output_directory  <- file.path(file_directory, "5) Output")
plots_output_directory <- file.path(file_output_directory, "Figures")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(R_objects_directory, "09) Process demultiplexed PacBio reads.RData"))
load(file.path(R_objects_directory, "10) Process demultiplexed reads - with subsampling.RData"))





# Set up loop -------------------------------------------------------------

DrawAllSmrtLinkPlots(DrawAccuracyHeatmap,
                     "Heatmaps",
                     "Accuracy heatmap"
                     )

DrawAllSubsampledPlots(DrawAccuracyHeatmap,
                       "Heatmaps",
                       "Accuracy heatmaps"
                       )




