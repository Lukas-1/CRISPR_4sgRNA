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
manuscript_directory     <- file.path(plots_output_directory, "Manuscript")



# Load data ---------------------------------------------------------------

load(file.path(s2rI_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2rI_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(s2rI_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))



# Prepare summary data frames ---------------------------------------------

use_summary_df <- ccs3_df_list[["original_summary_df"]]

CRISPRa_plate_numbers  <- plates_df[grepl("^HA", plates_df[, "Plate_name"]), "Plate_number"]
CRISPRko_plate_numbers <- plates_df[grepl("^HO", plates_df[, "Plate_name"]), "Plate_number"]

CRISPRa_summary_df  <- use_summary_df[use_summary_df[, "Plate_number"] %in% CRISPRa_plate_numbers, ]
CRISPRko_summary_df <- use_summary_df[use_summary_df[, "Plate_number"] %in% CRISPRko_plate_numbers, ]



# Draw example sand charts ------------------------------------------------

sg_sequences_df <- library_df[library_df[, "Modality"] %in% "CRISPRa", ]
DrawReorderedSandPlots(CRISPRa_summary_df)

sg_sequences_df <- library_df[library_df[, "Modality"] %in% "CRISPRko", ]
summary_df <- CRISPRko_summary_df[CRISPRko_summary_df[["Plate_number"]] %in% 109, ]
sg_sequences_df <- sg_sequences_df[sg_sequences_df[["Plate_number"]] %in% 109, ]
DrawReorderedSandPlots(summary_df)
rm(summary_df)

rm(sg_sequences_df)

SingleSandPlot(use_summary_df, invert_x_axis = TRUE)
SingleSandPlot(use_summary_df, invert_x_axis = TRUE, rotate_axes = FALSE, side_legend = FALSE)
SingleSandPlot(use_summary_df, invert_x_axis = TRUE, rotate_axes = TRUE, side_legend = FALSE)



# Draw example eCDF plots -------------------------------------------------

Plot_eCDF(CRISPRa_summary_df, "Count_at_least_1")
Plot_eCDF(CRISPRa_summary_df, "Count_all_4")
Plot_eCDF(CRISPRa_summary_df, "Count_all_4_promoters")

Plot_eCDF(CRISPRa_summary_df, "Num_reads_with_deletions_exceeding_20bp", flip_axis = TRUE)
Plot_eCDF(CRISPRa_summary_df, "4_guides", flip_axis = TRUE)



# Look up quantiles for cutoffs, or values for specific quantiles ---------

ValuesForQuantiles(CRISPRa_summary_df, "Count_at_least_1", 0.5)
ValuesForQuantiles(CRISPRa_summary_df, "Count_at_least_2", 0.5)
ValuesForQuantiles(CRISPRa_summary_df, "Count_at_least_3", 0.5)
ValuesForQuantiles(CRISPRa_summary_df, "Count_all_4",      0.5)

ValuesForQuantiles(CRISPRa_summary_df, "Count_at_least_1", 0.05)
ValuesForQuantiles(CRISPRa_summary_df, "Count_at_least_2", 0.05)
ValuesForQuantiles(CRISPRa_summary_df, "Count_at_least_3", 0.05)
ValuesForQuantiles(CRISPRa_summary_df, "Count_all_4",      0.05)


ValuesForQuantiles(CRISPRa_summary_df, "Count_at_least_3", 0.5)
ValuesForQuantiles(CRISPRa_summary_df, "Count_all_4", 0.5)
ValuesForQuantiles(CRISPRa_summary_df, "Count_at_least_3", 0.05)
ValuesForQuantiles(CRISPRa_summary_df, "Count_all_4", 0.05)
ValuesForQuantiles(CRISPRa_summary_df, "Num_reads_with_deletions_exceeding_20bp", 0.5)
FractionsForCutoffs(CRISPRa_summary_df, "Num_reads_with_deletions_exceeding_20bp", 0.5)

ValuesForQuantiles(CRISPRko_summary_df, "Count_at_least_3", 0.5)
ValuesForQuantiles(CRISPRko_summary_df, "Count_all_4", 0.5)
ValuesForQuantiles(CRISPRko_summary_df, "Count_at_least_3", 0.05)
ValuesForQuantiles(CRISPRko_summary_df, "Count_all_4", 0.05)


FractionsForCutoffs(CRISPRa_summary_df, "Count_all_4", 0.75)
FractionsForCutoffs(CRISPRa_summary_df, "All4_sg_cr_pass", 0.75)


median(ColumnToCDFVec(CRISPRa_summary_df, "Num_reads_with_deletions_spanning_tracrRNAs"))
median(ColumnToCDFVec(CRISPRa_summary_df, "Num_reads_with_deletions_spanning_promoters"))


table(c(CRISPRa_summary_df[, "Count_total"], CRISPRko_summary_df[, "Count_total"]) >= 10) /
(nrow(CRISPRa_summary_df) + nrow(CRISPRko_summary_df))


ValuesForQuantiles(CRISPRko_summary_df, "Count_sg1_cr1", 0.5)
ValuesForQuantiles(CRISPRko_summary_df, "Count_sg2_cr2", 0.5)
ValuesForQuantiles(CRISPRko_summary_df, "Count_sg3_cr3", 0.5)
ValuesForQuantiles(CRISPRko_summary_df, "Count_sg4_cr4", 0.5)

ValuesForQuantiles(CRISPRa_summary_df, "Count_sg1_cr1", 0.5)
ValuesForQuantiles(CRISPRa_summary_df, "Count_sg2_cr2", 0.5)
ValuesForQuantiles(CRISPRa_summary_df, "Count_sg3_cr3", 0.5)
ValuesForQuantiles(CRISPRa_summary_df, "Count_sg4_cr4", 0.5)



# Export sand charts for the manuscript -----------------------------------

manuscript_mai <- c(0.35 + 0.024, 0.5, 0.2 + 0.024, 0.9)

for (include_promoters in c(FALSE, TRUE)) {

  for (show_library in c("a", "o")) {

    if (show_library == "a") {
      current_summary_df <- CRISPRa_summary_df
      use_top_title <- "T.gonfio library"
    } else if (show_library == "o") {
      current_summary_df <- CRISPRko_summary_df
      use_top_title <- "T.spiezzo library"
    }
    figure_prefix <- "B"
    if (include_promoters) {
      if (show_library == "a") {
        figure_prefix <- "C"
      } else {
        figure_prefix <- "D"
      }
    }
    use_file_name <- paste0(figure_prefix, ") Sand chart - CRISPR", show_library, " library.pdf")
    pdf(file = file.path(manuscript_directory,
                         if (include_promoters) "Fig. S4" else "Fig. 4",
                         "Individual plots",
                         use_file_name
                         ),
        width = 3.4, height = 1.75
        )
    old_par <- par(cex = 0.7, lwd = 0.8, mai = manuscript_mai)
    SingleSandPlot(current_summary_df,
                   top_title             = use_top_title,
                   show_grid             = TRUE,
                   consider_promoters    = include_promoters,
                   rotate_axes           = FALSE,
                   invert_x_axis         = TRUE,
                   side_legend           = TRUE,
                   NA_in_legend          = FALSE,
                   set_mar               = FALSE,
                   data_axis_label       = "Percentage of reads",
                   vertical_y_label_line = 2.65,
                   x_label_line          = 2.1,
                   x_ticks_line          = 0.3,
                   use_mtext             = TRUE,
                   legend_x_start        = 1.3,
                   legend_x_title        = -0.1,
                   legend_y_gap          = 0.9
                   )
    par(old_par)
    dev.off()
  }
}



# Export eCDF plots for the manuscript ------------------------------------

for (include_promoters in c(FALSE, TRUE)) {

  for (show_library in c("a", "o")) {

    if (show_library == "a") {
      current_summary_df <- CRISPRa_summary_df
      use_top_title <- "T.gonfio library"
    } else if (show_library == "o") {
      current_summary_df <- CRISPRko_summary_df
      use_top_title <- "T.spiezzo library"
    }


    for (draw_figure in c("D", "E", "F")) {

      figure_prefix <- paste0(draw_figure, ") ")
      if (include_promoters) {
        if (draw_figure == "D") {
          figure_prefix <- file.path("Not used", figure_prefix)
        } else {
          next
        }
      }

      if (draw_figure == "D") {
        use_file_name <- paste0(figure_prefix, "eCDF - polyclonal bonus")
        use_line_colors <- c("#4E4073", brewer.pal(9, "Blues")[[6]])
        if (include_promoters) {
          column_combo <- "4_guides_with_promoters"
        } else {
          column_combo <- "4_guides"
        }
      } else if (draw_figure == "E") {
        use_file_name <- paste0(figure_prefix, "eCDF - deletions")
        column_combo <- "Deletions"
        use_line_colors <- brewer.pal(9, "Blues")[c(9, 7)]
        use_line_colors[[3]] <- colorRampPalette(c("#FDAC87", "#DB678B"))(10)[[4]] #brewer.pal(9, "Oranges")[[6]]
      } else if (draw_figure == "F") {
        use_file_name <- paste0(figure_prefix, "eCDF - contaminations")
        column_combo <- "Contaminations"
        use_line_colors <- brewer.pal(9, "Blues")[[7]] # brewer.pal(9, "Set1")[[7]]
      }

      use_file_name <- paste0(use_file_name, " - CRISPR", show_library, " library.pdf")
      pdf(file = file.path(manuscript_directory,
                           if (include_promoters) "Fig. S4" else "Fig. 4",
                           "Individual plots",
                           use_file_name
                           ),
          width = 3.4, height = 1.75
          )
      old_par <- par(cex = 0.7, lwd = 0.8, mai = manuscript_mai)
      Plot_eCDF(current_summary_df,
                column_combo,
                top_title             = use_top_title,
                flip_axis             = TRUE,
                rotate_axes           = FALSE,
                data_axis_label       = "Percentage of reads",
                set_mar               = FALSE,
                vertical_y_label_line = 2.65,
                x_label_line          = 2.1,
                x_ticks_line          = 0.3,
                use_mtext             = TRUE,
                always_side_legend    = TRUE,
                legend_y_gap          = if (draw_figure == "E") 1 else 1.125,
                legend_x_start        = if (draw_figure == "D") 0.6 else 0.7,
                legend_gap_ratio      = if (draw_figure == "E") 1.45 else 1.5,
                point_x_start         = 0.9,
                reverse_legend_order  = draw_figure == "D",
                line_colors           = use_line_colors,
                legend_pch            = 22,
                lwd_multiplier        = 1.5
                )
      par(old_par)
      dev.off()
    }
  }
}



# Export all eCDF plots ---------------------------------------------------

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




# Export individual graphics ----------------------------------------------

use_plate_numbers <- unique(library_df[["Plate_number"]])

DrawBarplotsAndHeatmapsForAllPlates(export_PNGs = FALSE)





