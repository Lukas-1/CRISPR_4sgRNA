### 28 May 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "09) Producing heatmaps.R"))
source(file.path(R_functions_directory, "11) Creating stacked barplots for visualizing alterations.R"))



# Define folder paths -----------------------------------------------------

sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")
file_output_directory    <- file.path(sql2_directory, "5) Output")
plots_output_directory   <- file.path(file_output_directory, "Figures")



# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(sql2_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(sql2_R_objects_directory, "09) Process demultiplexed PacBio reads.RData"))


# # try ---------------------------------------------------------------------
#
# use_df <- ccs7_df_list[["individual_reads_df"]]
# pass_filters <- use_df[["Passes_filters"]] == 1
#
# use_df <- use_df[pass_filters, ]
#
#
# use_df <- use_df[use_df[["Plate_number"]] %in% 13, ]
#
#
# sub_df <- use_df[use_df[["Well_number"]] %in% 2, ]
# table(sub_df[["Contam_well"]])
#
# sub_df <- use_df[use_df[["Well_number"]] %in% 4, ]
# table(sub_df[["Contam_well"]])
#
# sub_df <- use_df[use_df[["Well_number"]] %in% 6, ]
# table(sub_df[["Contam_well"]])
#
# sub_df <- use_df[use_df[["Well_number"]] %in% 8, ]
# table(sub_df[["Contam_well"]])
#
# sub_df <- use_df[use_df[["Well_number"]] %in% 10, ]
# table(sub_df[["Contam_well"]])
#
# sub_df <- use_df[use_df[["Well_number"]] %in% 25, ]
# table(sub_df[["Contam_well"]])
#
# sub_df <- use_df[use_df[["Well_number"]] %in% 26, ]
#
#
#
# vapply(seq(2, 382, by = 2), function(x) {
#   are_this_well <- use_df[["Well_number"]] %in% x
#   names(sort(table(use_df[["Contam_well"]][are_this_well]), decreasing = TRUE))[[1]]
# }, "")
#
#
# library_mat <- as.matrix(library_df[, paste0("Sequence_sg", 1:4)])
#
# "ACTCTGCGTTGAGAGAAACA" %in% library_mat
#
#
# use_df[use_df[["ZMW"]] %in% 80020138, ]
#
#
# groo_df <- extracted_df[extracted_df[["ZMW"]] %in% 80020138, ]
#
# groo_df[groo_df[["Feature"]] %in% paste0("sg", 1:4), ]
#



# Export individual graphics ----------------------------------------------

ccs_numbers <- c(3, 5, 7)
accuracy_percentages <- c(99, 99.9, 99.99)

use_plate_numbers <- plates_df[["Plate_number"]][order(plates_df[["Plate_rank"]])]


for (reorder_wells in c(TRUE, FALSE)) {

  order_folder <- c("original order", "re-ordered")[[as.integer(reorder_wells) + 1]]
  heatmaps_folder   <- paste0("Heatmaps - ", order_folder)
  barplots_folder   <- paste0("Stacked barplots - ", order_folder)
  sand_charts_folder <- "Sand charts"

  for (i in seq_along(ccs_numbers)) {

    use_df_list <- get(paste0("ccs", ccs_numbers[[i]], "_df_list"))

    titles_list <- lapply(plates_df[["Plate_number"]], function(x) {
      bquote({"Plate #" * .(as.character(x)) * " \u2013 " *
          bold(.(plates_df[x, "Plate_name"])) *
          " (" >= .(as.character(ccs_numbers[[i]])) * " consensus reads "} *
            "and " >= .(as.character(accuracy_percentages[[i]])) * "% accuracy)"
      )
    })

    for (filter_stage in 1:2) {

      df_name <- c("original_summary_df", "filtered_summary_df")[[filter_stage]] # "filtered_gRNAs_df"

      use_summary_df <- use_df_list[[df_name]]

      file_name <- paste0("CCS", ccs_numbers[[i]],
                          " (", accuracy_percentages[[i]], ") - ",
                          c("i) unfiltered", "ii) filtered", "iii) filtered gRNAs")[[filter_stage]]
                          )

      pdf(file = file.path(plots_output_directory, heatmaps_folder, paste0("Heatmaps - ", file_name, ".pdf")),
          width  = use_width, height = use_height
          )
      for (plate_number in use_plate_numbers) {
        sub_df <- use_summary_df[use_summary_df[["Plate_number"]] %in% plate_number, ]
        sg_sequences_df <- library_df[library_df[["Plate_number"]] %in% plate_number, ]
        DrawAccuracyHeatmap(sub_df,
                            main_title = titles_list[[plate_number]],
                            reorder_wells = reorder_wells
                            )
      }
      dev.off()

      if (reorder_wells) {
        pdf(file = file.path(plots_output_directory, sand_charts_folder, paste0("Sand chart - ", file_name, ".pdf")),
            width = 3.8, height = 6.7
            )
        for (plate_number in use_plate_numbers) {
          sub_df <- use_summary_df[use_summary_df[["Plate_number"]] %in% plate_number, ]
          sg_sequences_df <- library_df[library_df[["Plate_number"]] %in% plate_number, ]
          DrawReorderedSandPlots(sub_df,
                                 main_title = titles_list[[plate_number]]
                                 )
        }
        dev.off()
      }


      pdf(file = file.path(plots_output_directory, barplots_folder, paste0("Stacked barplots - ", file_name, ".pdf")),
          height = use_height,
          width  = use_width
          )
      for (plate_number in use_plate_numbers) {
        sub_df <- use_summary_df[use_summary_df[["Plate_number"]] %in% plate_number, ]
        sg_sequences_df <- library_df[library_df[["Plate_number"]] %in% plate_number, ]
        DrawAlterationBarplot(sub_df,
                              main_title = titles_list[[plate_number]],
                              reorder_wells = reorder_wells
                              )
      }
      dev.off()
    }
  }
}






