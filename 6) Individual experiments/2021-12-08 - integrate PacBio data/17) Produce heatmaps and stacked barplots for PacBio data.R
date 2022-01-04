### 31st December 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "09) Producing heatmaps.R"))
source(file.path(R_functions_directory, "11) Creating stacked barplots for visualizing alterations.R"))




# Define folder paths -----------------------------------------------------

s2rI_directory           <- file.path(experiments_directory, "2021-12-08 - integrate PacBio data")
s2rI_R_objects_directory <- file.path(s2rI_directory, "3) R objects")
file_output_directory    <- file.path(s2rI_directory, "5) Output")
plots_output_directory   <- file.path(file_output_directory, "Figures")
PNGs_output_directory    <- file.path(file_output_directory, "PNGs")
across_plate_directory   <- file.path(PNGs_output_directory, "Sand charts", "Summaries across plates")




# Load data ---------------------------------------------------------------

load(file.path(s2rI_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2rI_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(s2rI_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))





# Try stuff ---------------------------------------------------------------

use_summary_df <- ccs3_df_list[["original_summary_df"]]

matches_vec <- match(use_summary_df[, "Plate_number"], plates_df[, "Plate_number"])
plate_names_vec <- plates_df[matches_vec, "Plate_name"]
are_CRISPRa <- grepl("^HA", plate_names_vec)
are_CRISPRko <- grepl("^HO", plate_names_vec)

CRISPRa_summary_df <- use_summary_df[are_CRISPRa, ]
CRISPRko_summary_df <- use_summary_df[are_CRISPRko, ]
row.names(CRISPRa_summary_df) <- NULL
row.names(CRISPRko_summary_df) <- NULL


summary_df <- CRISPRa_summary_df[CRISPRa_summary_df[["Plate_number"]] %in% 47, ]

SingleSandPlot(summary_df)

for (show_grid in c(TRUE, FALSE)) {
  for (show_library in c("a", "o")) {
    for (increasing_order in c(FALSE, TRUE)) {

      file_name <- paste0("Sand chart - CRISPR", show_library, " library")
      if (increasing_order) {
        file_name <- paste0(file_name, " - increasing order")
        folder_path <- file.path(across_plate_directory, "Increasing order")
      } else {
        file_name <- paste0(file_name, " - decreasing order")
        folder_path <- file.path(across_plate_directory, "Decreasing order")
      }
      if (show_grid) {
        file_name <- paste0(file_name, " - with grid")
      }
      file_name <- paste0(file_name, ".png")
      if (show_library == "a") {
        current_summary_df <- CRISPRa_summary_df
      } else if (show_library == "o") {
        current_summary_df <- CRISPRko_summary_df
      }
      x_axis_label <- paste0("% wells in the CRISPR", show_library, " library")
      png(file = file.path(folder_path, file_name),
          width = 6.3, height = 4.5, units = "in", res = 600
          )
      SingleSandPlot(current_summary_df,
                     show_grid = show_grid,
                     invert_x_axis = !(increasing_order),
                     x_axis_label = x_axis_label
                     )
      dev.off()

    }
  }
}



sg_sequences_df <- library_df[library_df[, "Modality"] %in% "CRISPRa", ]
DrawReorderedSandPlots(CRISPRa_summary_df)

PlotBarplotMat(column_mat_list[[4]], use_colors, positions_vec)





# Export individual graphics ----------------------------------------------

use_plate_numbers <- unique(library_df[["Plate_number"]])

DrawBarplotsAndHeatmapsForAllPlates(export_PNGs = FALSE)







