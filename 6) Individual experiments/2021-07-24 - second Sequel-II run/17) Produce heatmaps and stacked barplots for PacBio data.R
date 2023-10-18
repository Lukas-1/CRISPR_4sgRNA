### 29th July 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "09) Producing heatmaps.R"))
source(file.path(R_functions_directory, "11) Creating stacked barplots for visualizing alterations.R"))




# Define folder paths -----------------------------------------------------

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
file_output_directory    <- file.path(s2r2_directory, "5) Output")
plots_output_directory   <- file.path(file_output_directory, "Figures")
PNGs_output_directory    <- file.path(file_output_directory, "PNGs")




# Load data ---------------------------------------------------------------

load(file.path(s2r2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2r2_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(s2r2_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))





# Export individual graphics ----------------------------------------------

use_plate_numbers <- plates_df[["Plate_number"]]

DrawBarplotsAndHeatmapsForAllPlates()




# Export .emf files from a selected plate ---------------------------------

library("devEMF")
thesis_dir <- file.path(plots_output_directory, "Thesis")


summary_df <- ccs7_df_list[["filtered_summary_df"]]
summary_df <- summary_df[summary_df[, "Plate_number"] %in% 33, ]
row.names(summary_df) <- NULL


use_width <- 3.0
use_height <- 1.5
use_lwd <- 0.9
use_cex <- 0.7

emf(file = file.path(thesis_dir, paste0("Plate HA_21 - stacked barplot.emf")),
    width = use_width, height = use_height,
    emfPlus = FALSE
    )
par(mai = c(0.4, 0.4, 0.13, 0.13), lwd = use_lwd, cex = use_cex)
MakeEmptyPlot()
DrawPercentCorrectBarplot(summary_df,
                          prefix_text      = "",
                          postfix_text     = " correct gRNAs",
                          text_at_bottom   = TRUE,
                          narrow_lwd       = FALSE,
                          color_box_legend = TRUE,
                          aspect_ratio     = (use_width  - sum(par("mai")[c(2, 4)])) /
                                              (use_height - sum(par("mai")[c(1, 3)]))
                          )

dev.off()


emf(file = file.path(thesis_dir, paste0("Plate HA_21 - heatmap.emf")),
    width = use_width, height = use_height,
    emfPlus = FALSE
    )
par(mai = c(0.4, 0.4, 0.13, 0.13), lwd = use_lwd, cex = use_cex)
DrawHeatMap(summary_df,
            trapezoid_on_new_plot = FALSE,
            trapezoid_start_y     = -0.16,
            trapezoid_y           = -0.125,
            trapezoid_end_y       = -0.05,
            trapezoid_start_x     = 0.75,
            trapezoid_end_x       = 1,
            add_percent           = FALSE,
            accuracy_label        = "% correct",
            accuracy_label_cex    = 1,
            label_y_factor        = 1.15,
            use_lwd               = 0.75,
            bold_percentages      = FALSE
            )
dev.off()





