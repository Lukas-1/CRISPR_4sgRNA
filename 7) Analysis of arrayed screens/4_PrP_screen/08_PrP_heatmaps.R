# 2022-01-20


# Load packages and source code -------------------------------------------

library("squash")
library("viridis")
library("RColorBrewer")

project_dir   <- "~/R_projects/CRISPRa_TF"
functions_dir <- file.path(project_dir, "1_R_scripts", "1_R_functions")
source(file.path(functions_dir, "1_General_functions", "01_labels_and_annotations.R"))
source(file.path(functions_dir, "1_General_functions", "02_plotting_helper_functions.R"))
source(file.path(functions_dir, "2_Analyzing_data",    "01_calculating_scores.R"))
source(file.path(functions_dir, "3_Visualizing_data",  "06_Heatmaps.R"))



# Define folder path ------------------------------------------------------

r_data_dir <- file.path(project_dir, "3_R_objects", "3_PrP")
output_dir <- file.path(project_dir, "4_output", "PrP")



# Load data ---------------------------------------------------------------

load(file.path(r_data_dir, "02_analyse_data.RData"))



# Modify labels (for PrPc screen) -----------------------------------------

AdjustLabels()



# Draw example heatmaps ---------------------------------------------------

HeatmapForPlate(PrP_df, 10, "Raw_rep1")
HeatmapForPlate(PrP_df, 10, "Raw_rep2")


HeatmapForPlate(PrP_df, 10, "Raw_rep1", weighting_for_controls = FALSE,
                use_one_scale = FALSE
                )
HeatmapForPlate(PrP_df, 10, "Raw_rep1", weighting_for_controls = FALSE,
                use_one_scale = TRUE
                )

HeatmapForPlate(PrP_df, 12, "p_value_log2",
                use_subtext = long_column_labels[["p_value_log2"]]
                )

AveragedHeatmap(PrP_df, "Raw_rep1", both_replicates = TRUE,  use_one_scale = TRUE)
AveragedHeatmap(PrP_df, "Raw_rep1", both_replicates = TRUE,  use_one_scale = FALSE)
AveragedHeatmap(PrP_df, "Raw_rep1", both_replicates = FALSE, use_one_scale = TRUE)
AveragedHeatmap(PrP_df, "Raw_rep1", both_replicates = FALSE, use_one_scale = FALSE)
AveragedHeatmap(PrP_df, "Raw_rep2", both_replicates = FALSE, use_one_scale = TRUE)
AveragedHeatmap(PrP_df, "Raw_rep2", both_replicates = FALSE, use_one_scale = FALSE)
AveragedHeatmap(PrP_df, "Raw_rep2", both_replicates = FALSE, use_one_scale = TRUE, weighting_for_controls = FALSE)

AveragedHeatmap(PrP_df, "Raw_rep1", both_replicates = TRUE,  use_one_scale = TRUE, take_median = TRUE)



HeatmapForPlate(PrP_df, 1, "Raw_rep1", use_one_scale = TRUE,  weighting_for_controls = TRUE) # default
HeatmapForPlate(PrP_df, 1, "Raw_rep1", use_one_scale = FALSE, weighting_for_controls = TRUE)
HeatmapForPlate(PrP_df, 1, "Raw_rep1", use_one_scale = TRUE,  weighting_for_controls = FALSE)
HeatmapForPlate(PrP_df, 1, "Raw_rep1", use_one_scale = FALSE, weighting_for_controls = FALSE)



HeatmapForPlate(PrP_df, 1, "Raw_rep1",
                ColorFunction = magma
                )


HeatmapForPlate(PrP_df, 1, "Raw_log2_rep1",
                weighting_for_controls = FALSE,
                use_one_scale = FALSE,
                use_subtext = long_column_labels[["Raw_log2_rep1"]]
                )



HeatmapForPlate(PrP_df, 2, "CellTiterGlo_raw")
HeatmapForPlate(PrP_df, 2, "CellTiterGlo_raw", uniform_legend = FALSE)


HeatmapForPlate(PrP_df, 2, "CellTiterGlo_foldNT")
HeatmapForPlate(PrP_df, 2, "CellTiterGlo_foldNT",
                weighting_for_controls = TRUE
                )

HeatmapForPlate(PrP_df, 2, "CellTiterGlo_foldNT",
                num_uniform_breaks = 100
                )

HeatmapForPlate(PrP_df, 12, "CellTiterGlo_raw",
                weighting_for_controls = TRUE
                )
HeatmapForPlate(PrP_df, 12, "CellTiterGlo_raw",
                weighting_for_controls = TRUE,
                num_other_breaks = 200
                )

HeatmapForPlate(PrP_df, 7, "CellTiterGlo_raw",
                weighting_for_controls = TRUE
                )

HeatmapForPlate(PrP_df, 5, "CellTiterGlo_raw",
                weighting_for_controls = TRUE
                )


HeatmapForPlate(PrP_df, 1, "SSMD_deltaNT")
HeatmapForPlate(PrP_df, 1, "SSMD_log2")
HeatmapForPlate(PrP_df, 1, "SSMD_deltaNT_Glo", uniform_legend = TRUE)
HeatmapForPlate(PrP_df, 6, "SSMD_deltaNT")


HeatmapForPlate(PrP_df, 10, "SSMD_log2")
HeatmapForPlate(PrP_df, 2, "p_value_log2", num_uniform_breaks = 100,
                use_subtext = expression(italic("p") * " value")
                )

HeatmapForPlate(PrP_df, 1, "Raw_rep1")
HeatmapForPlate(PrP_df, 2, "Raw_rep2")


HeatmapForPlate(PrP_df, 1, "DeltaNT_rep1")
HeatmapForPlate(PrP_df, 1, "DeltaNT_rep2")


HeatmapForPlate(PrP_df, 10, "Raw_rep1")
HeatmapForPlate(PrP_df, 10, "Raw_rep2")


HeatmapForPlate(PrP_df, 1, "SSMD_log2", label_values = TRUE)
HeatmapForPlate(PrP_df, 1, "Raw_rep1", label_values = TRUE)


HeatmapForPlate(PrP_df, 1, "FoldNT_rep1")
HeatmapForPlate(PrP_df, 1, "FoldNT_rep1", take_log2 = TRUE)


HeatmapForPlate(PrP_df, 5, "Hit_strength_deltaNT")


HeatmapForPlate(PrP_df, 10, "Raw_rep1", uniform_legend = TRUE)
HeatmapForPlate(PrP_df, 10, "Raw_rep1", uniform_legend = FALSE)

HeatmapForPlate(PrP_df, 10, "Raw_rep2")


HeatmapForPlate(PrP_df, 1, "SSMD_log2", label_values = TRUE)
HeatmapForPlate(PrP_df, 1, "Raw_rep1", label_values = TRUE)


HeatmapForPlate(PrP_df, 1, "FoldNT_rep1")
HeatmapForPlate(PrP_df, 1, "FoldNT_rep1", take_log2 = TRUE)

HeatmapForPlate(PrP_df, 1, "Raw_log2_rep1")
HeatmapForPlate(PrP_df, 1, "Raw_log2_rep1", weighting_for_controls = FALSE)




# Export heatmaps as PDF or PNG files -------------------------------------

if (TRUE) {
  ExportAllHeatmaps(PrP_df)
}




