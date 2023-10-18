### 28 May 2020 ###




# Import packages and source code -----------------------------------------

CRISPR_root_directory  <- "~/CRISPR_4sgRNA"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "09) Producing heatmaps.R"))




# Define folder paths -----------------------------------------------------

R_objects_directory    <- file.path(file_directory, "3) R objects")
file_output_directory  <- file.path(file_directory, "5) Output")
plots_output_directory <- file.path(file_output_directory, "Figures")
manuscript_directory   <- file.path(plots_output_directory, "For the manuscript")



# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(R_objects_directory, "11) Process demultiplexed PacBio reads.RData"))
load(file.path(R_objects_directory, "12) Process demultiplexed reads - with subsampling.RData"))




# Export individual graphics ----------------------------------------------

use_summary_df <- PrepareSummaryDf(sl7_ccs7_df_list[["filtered_summary_df"]])

ExportFiguresForManuscript(use_summary_df, "First 384-well plate")




# Set up loop -------------------------------------------------------------

AllSmrtLinkVersionsForOnePlate(DrawAccuracyHeatmap,
                               "Heatmaps",
                               "Accuracy heatmap"
                               )

DrawAllSubsampledPlotsForOnePlate(DrawAccuracyHeatmap,
                                  "Heatmaps",
                                  "Accuracy heatmaps"
                                  )




