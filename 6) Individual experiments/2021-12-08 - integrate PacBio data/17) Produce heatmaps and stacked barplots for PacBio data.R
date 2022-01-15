### 31st December 2021 ###



# Import packages and source code -----------------------------------------

library("matrixStats")
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
# PNGs_output_directory    <- file.path(file_output_directory, "PNGs")
across_plate_directory   <- file.path(plots_output_directory, "Sand and eCDF charts - all plates combined")




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


Plot_eCDF(CRISPRa_summary_df, "Count_at_least_1")

Plot_eCDF(CRISPRa_summary_df, "Count_all_4")
Plot_eCDF(CRISPRa_summary_df, "Count_all_4_promoters")



summary_df <- CRISPRa_summary_df[CRISPRa_summary_df[["Plate_number"]] %in% 47, ]


SingleSandPlot(summary_df, invert_x_axis = FALSE, rotate_axes = TRUE)


try_steps <- MakeSteps(delete_column_mat[, 3])
plot(y = try_steps[, "data_axis"], x = try_steps[, "fraction_axis"], xlim = c(0, 1), ylim = c(0, 1), type = "l")
DrawPolygon(AddCorner(try_steps), rotate_axes = TRUE)

pdf(file = "foo.pdf", width = 5, height = 5)
SingleSandPlot(summary_df, invert_x_axis = FALSE)
dev.off()

# sliced_mat <- SliceFraction(try_steps, "fraction_axis", 0.98, 0.995)
# DrawPolygon(sliced_mat, fill_color = "red", rotate_axes = FALSE, flip_axis = TRUE)


basic_rectangle <- cbind("data_axis"     = c(0, 1, 1),
                         "fraction_axis" = c(0, 0, 1)
                         )


sliced_mat <- SliceFraction(basic_rectangle, "fraction_axis", 0.2, 0.4)
DrawPolygon(sliced_mat, "red")

sliced_mat <- SliceFraction(try_steps, "fraction_axis", 0, 0.02)
DrawPolygon(sliced_mat, fill_color = "red", rotate_axes = FALSE, flip_axis = FALSE)


sliced_mat <- SliceFraction(try_steps, "data_axis", 0.05, 0.1)
DrawPolygon(sliced_mat, fill_color = "blue", rotate_axes = FALSE, flip_axis = FALSE)


for (i in seq_len(nrow(sliced_mat))) {
  points(x = sliced_mat[i, 1], y = sliced_mat[i, 2], pch = 16, xpd = NA, col = "blue", cex = 1)
}



are_same <- round(try_steps[, "fraction_axis"], digits = 12) == sort(round(try_steps[, "fraction_axis"], digits = 12))


try_steps <- MakeSteps(delete_column_mat[, 1])
plot(y = try_steps[, "data_axis"], x = try_steps[, "fraction_axis"], xlim = c(0, 1), ylim = c(0, 1), type = "l")
polygon(y = try_steps[, "data_axis"], x = try_steps[, "fraction_axis"], xlim = c(0, 1), ylim = c(0, 1), col = "blue")


SingleSandPlot(summary_df, invert_x_axis = TRUE)

MakeEmptyPlot()
DrawPolygon(try_steps, rotate_axes = TRUE)

MakeEmptyPlot()
DrawPolygon(try_steps, rotate_axes = FALSE, flip_axis = TRUE)

polygon(y = try_steps[, "data_axis"], x = try_steps[, "fraction_axis"])



rotate_axes_vec      <- c(FALSE, FALSE, TRUE, TRUE)
increasing_order_vec <- c(FALSE, TRUE, FALSE, TRUE)
folder_names_vec     <- paste0(letters[1:4], ") Percent correct on ",
                               ifelse(rotate_axes_vec, "y", "x"), " axis - ",
                               ifelse(increasing_order_vec, "increasing", "decreasing"),
                               " order"
                               )


for (i in 1:4) {

  rotate_axes <- rotate_axes_vec[[i]]
  increasing_order <- increasing_order_vec[[i]]
  folder_path <- file.path(across_plate_directory, folder_names_vec[[i]])

  for (include_promoters in c(FALSE, TRUE)) {

    if (include_promoters) {
      promoter_prefix <- " - with promoters"
    } else {
      promoter_prefix <- " - only sg_cr"
    }

   for (show_library in c("a", "o")) {

      file_suffix <- paste0(" - CRISPR", show_library, " library.pdf")
      if (show_library == "a") {
        current_summary_df <- CRISPRa_summary_df
      } else if (show_library == "o") {
        current_summary_df <- CRISPRko_summary_df
      }
      fraction_axis_label <- paste0("Percentile of wells in the CRISPR", show_library, " library")
      pdf(file = file.path(folder_path, paste0("Sand chart", promoter_prefix, file_suffix)),
          width = 6.3, height = 4.5
          )
      SingleSandPlot(current_summary_df,
                     show_grid = TRUE,
                     consider_promoters = include_promoters,
                     rotate_axes = rotate_axes,
                     invert_x_axis = !(increasing_order),
                     fraction_axis_label = fraction_axis_label
                     )
      SingleSandPlot(current_summary_df,
                     show_grid = FALSE,
                     consider_promoters = include_promoters,
                     rotate_axes = rotate_axes,
                     invert_x_axis = !(increasing_order),
                     fraction_axis_label = fraction_axis_label
                     )
      dev.off()

      if (include_promoters) {
        pdf(file = file.path(folder_path, paste0("eCDF plots", file_suffix)),
            width = 7.2, height = 5.5
            )
        for (use_column in percentages_metrics) {
          Plot_eCDF(current_summary_df,
                    use_column,
                    rotate_axes         = rotate_axes,
                    flip_axis           = !(increasing_order),
                    fraction_axis_label = fraction_axis_label
                    )
        }
        dev.off()
      }
    }
  }
}




sg_sequences_df <- library_df[library_df[, "Modality"] %in% "CRISPRa", ]
DrawReorderedSandPlots(CRISPRa_summary_df)



sg_sequences_df <- library_df[library_df[, "Modality"] %in% "CRISPRko", ]
summary_df <- CRISPRko_summary_df[CRISPRko_summary_df[["Plate_number"]] %in% 109, ]
sg_sequences_df <- sg_sequences_df[sg_sequences_df[["Plate_number"]] %in% 109, ]
DrawReorderedSandPlots(summary_df)






# Export individual graphics ----------------------------------------------

use_plate_numbers <- unique(library_df[["Plate_number"]])

DrawBarplotsAndHeatmapsForAllPlates(export_PNGs = FALSE)







