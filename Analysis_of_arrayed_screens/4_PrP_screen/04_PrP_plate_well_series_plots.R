# 2022-01-18


# Load packages and source code -------------------------------------------

library("RColorBrewer")

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_scripts", "1_R_functions")
source(file.path(functions_dir, "1_General_functions", "01_labels_and_annotations.R"))
source(file.path(functions_dir, "1_General_functions", "02_plotting_helper_functions.R"))
source(file.path(functions_dir, "3_Visualizing_data",  "02_Plate_well_series_plots.R"))



# Define folder path ------------------------------------------------------

r_data_dir <- file.path(project_dir, "3_R_objects", "3_PrP")
output_dir <- file.path(project_dir, "4_output", "PrP")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "02_analyse_data.RData"))



# Modify labels (for PrPc screen) -----------------------------------------

AdjustLabels()
controls_labels[["Pos"]] <- c("Positive", "controls", expression("(" * italic("PRNP") * " gene)"))




# Draw example plots ------------------------------------------------------

PlateWellPlot(PrP_df, emphasize_NT = TRUE)
PlateWellPlot(PrP_df, order_by_column = TRUE)

PlateWellPlot(PrP_df, order_by_column = TRUE, aggregate_wells = TRUE)
PlateWellPlot(PrP_df, order_by_column = FALSE, aggregate_wells = FALSE)


PlateWellPlot(PrP_df, "Raw_Glo_rep1", order_by_column = TRUE, aggregate_wells = TRUE)
PlateWellPlot(PrP_df, "Raw_Glo_rep1", order_by_column = FALSE, aggregate_wells = TRUE)

PlateWellPlot(PrP_df, "FoldNT_rep1")
PlateWellPlot(PrP_df, "CellTiterGlo_raw")

PlateWellPlot(PrP_df, "DeltaNT_rep1")
PlateWellPlot(PrP_df, "Raw_log2_rep1")

PlateWellPlot(PrP_df, "Hit_strength_deltaNT_Glo")




# Export plots as PDF and PNG files ---------------------------------------

series_top_folder <- file.path(output_dir, "Figures", "Plate well series plots")

ExportAllPlateSeriesPlots(PrP_df, top_folder = series_top_folder)



