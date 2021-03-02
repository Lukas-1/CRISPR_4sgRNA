### 21st October 2020 ###




# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "09) Producing heatmaps.R"))




# Define folder paths -----------------------------------------------------

plate2_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-18 - PacBio - second 384-well plate")
p2_R_objects_directory <- file.path(plate2_directory, "2) R objects")
file_output_directory  <- file.path(plate2_directory, "3) Output")
plots_output_directory <- file.path(file_output_directory, "Figures")
manuscript_directory   <- file.path(plots_output_directory, "For the manuscript")




# Load data ---------------------------------------------------------------

load(file.path(p2_R_objects_directory, "01) Import and process sgRNA sequences.RData"))
load(file.path(p2_R_objects_directory, "07) Process demultiplexed PacBio reads.RData"))
load(file.path(p2_R_objects_directory, "08) Process demultiplexed reads - with subsampling.RData"))






# Exclude problematic wells, for the time being ---------------------------

sg_sequences_df[["Empty_well"]] <- ifelse(sg_sequences_df[["Well_number"]] == 2,
                                          TRUE, FALSE
                                          )





# Export individual graphics ----------------------------------------------

use_summary_df <- PrepareSummaryDf(sl7_ccs7_df_list[["filtered_summary_df"]])
use_summary_df <- use_summary_df[use_summary_df[["Block"]] %in% 1, ]

ExportFiguresForManuscript(use_summary_df, "Colony-picked controls")





# Set up loop -------------------------------------------------------------

DrawAllSmrtLinkPlots(DrawAccuracyHeatmap,
                     "Heatmaps",
                     "Accuracy heatmap"
                     )

DrawAllSubsampledPlots(DrawAccuracyHeatmap,
                       "Heatmaps",
                       "Accuracy heatmaps"
                       )





