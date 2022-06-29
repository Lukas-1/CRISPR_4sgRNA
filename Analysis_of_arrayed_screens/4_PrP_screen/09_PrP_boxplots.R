# 2022-02-11


# Load packages and source code -------------------------------------------

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_scripts", "1_R_functions")
source(file.path(functions_dir, "1_General_functions", "01_labels_and_annotations.R"))
source(file.path(functions_dir, "1_General_functions", "02_plotting_helper_functions.R"))
source(file.path(functions_dir, "3_Visualizing_data",  "07_Box_plots.R"))



# Define folder path ------------------------------------------------------

r_data_dir  <- file.path(project_dir, "3_R_objects", "3_PrP")
output_dir  <- file.path(project_dir, "4_output", "PrP")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "02_analyse_data.RData"))




# Modify labels (for PrPc screen) -----------------------------------------

AdjustLabels()



# Draw example plots ------------------------------------------------------

BeeBoxPlates(PrP_df, "Raw_rep1", split_NT = TRUE)

BeeBoxPlates(PrP_df, "Raw_rep1", show_subgroups = TRUE, plate_number = 2)
BeeBoxPlates(PrP_df, "Raw_rep1", show_subgroups = TRUE, show_96wp = TRUE)
BeeBoxPlates(PrP_df, "CellTiterGlo_raw", show_subgroups = TRUE, show_96wp = TRUE, plate_number = 22)
BeeBoxPlates(PrP_df, "CellTiterGlo_raw", show_subgroups = FALSE, show_96wp = TRUE, plate_number = 22)

BeeBoxPlates(PrP_df, "Raw_Glo_rep1")

BeeBoxPlates(PrP_df, "Raw_rep1", compare_group = "Gene")
BeeBoxPlates(PrP_df, "Raw_rep1", compare_group = "Own NT control")
BeeBoxPlates(PrP_df, "Raw_rep1", compare_group = "Pos. control")
BeeBoxPlates(PrP_df, "Raw_rep1", compare_group = "Genes and NT")

BeeBoxPlates(PrP_df, "Raw_rep1", compare_group = "Gene", show_96wp = TRUE)
BeeBoxPlates(PrP_df, "Raw_rep1", compare_group = "NT control", show_96wp = TRUE)





# Export plots as PDF and PNG files ---------------------------------------

ExportAllBoxPlots(PrP_df, file.path(output_dir, "Figures", "Box plots"))






