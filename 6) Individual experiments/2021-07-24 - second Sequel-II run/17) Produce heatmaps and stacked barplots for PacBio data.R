### 29th July 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
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
# load(file.path(s2r2_R_objects_directory, "13) Process demultiplexed reads - with subsampling.RData"))




# Export individual graphics ----------------------------------------------

use_plate_numbers <- plates_df[["Plate_number"]]

DrawBarplotsAndHeatmapsForAllPlates(export_PNGs = FALSE)






# # Export the results of subsampling simulations ---------------------------
#
# sampling_levels <- names(subsampled_list)
# are_all <- sampling_levels == "100% sampled"
# title_prefixes <- ifelse(are_all,
#                          "All reads",
#                          paste0("Sampled ", sub(" sampled", "", sampling_levels, fixed = TRUE), " of reads")
#                          )
# use_num_reps <- 3L
#
# for (file_format in c("png", "pdf")) {
#   for (i in seq_along(ccs_numbers)) {
#     for (filter_stage in 1:2) {
#
#       df_name <- c("original_summary_df", "filtered_summary_df")[[filter_stage]] # "filtered_gRNAs_df"
#       use_summary_df <- use_df_list[[df_name]]
#       sel_name <- paste0("CCS", ccs_numbers[[i]],
#                          " (", accuracy_percentages[[i]], ") - ",
#                          c("i) unfiltered", "ii) filtered", "iii) filtered gRNAs")[[filter_stage]]
#                          )
#
#       for (make_plot in c("heatmap", "sparse heatmap", "stacked barplot")) {
#
#         if (make_plot == "heatmap") {
#           folder_name <- "Heatmaps"
#           PlotFunction <- DrawAccuracyHeatmap
#           current_width <- use_width
#           current_height <- use_height
#         } else if (make_plot == "sparse heatmap") {
#           folder_name <- "Heatmaps (sparse)"
#           PlotFunction <- function(summary_df, main_title, reorder_wells) SparseAccuracyHeatmap(summary_df, main_title)
#           current_width <- sparse_width
#           current_height <- sparse_height
#         } else if (make_plot == "stacked barplot") {
#           folder_name <- "Stacked barplots"
#           PlotFunction <- ModifiedAlterationBarplot
#           current_width <- alterations_width
#           current_height <- alterations_height
#         }
#
#         folder_path <- file.path("Sub-sampled", folder_name, sel_name)
#         if (file_format == "pdf") {
#           folder_path <- file.path(plots_output_directory, folder_path)
#         } else if (file_format == "png") {
#           folder_path <- file.path(PNGs_output_directory, folder_path)
#         }
#         dir.create(folder_path, showWarnings = FALSE)
#
#         for (plate_number in use_plate_numbers) {
#           plate_name <- paste0("Plate", formatC(plate_number, width = 2, flag = "0"),
#                                " - ",
#                                plates_df[["Plate_name"]][plates_df[["Plate_number"]] == plate_number]
#                                )
#
#           sg_sequences_df <- library_df[library_df[["Plate_number"]] %in% plate_number, ]
#
#           if (file_format == "pdf") {
#             pdf(file = file.path(folder_path, paste0(plate_name, ".pdf")),
#                 width = current_width, height = current_height
#                 )
#           } else if (file_format == "png") {
#             sub_folder_path <- file.path(folder_path, plate_name)
#             dir.create(sub_folder_path, showWarnings = FALSE)
#           }
#           for (j in seq_along(subsampled_list)) {
#             for (rep in seq_len(min(length(subsampled_list[[j]]), use_num_reps))) {
#               title_prefix <- title_prefixes[[j]]
#               is_all <- names(subsampled_list)[[j]] == "100% sampled"
#               if (!(is_all)) {
#                 title_prefix <- paste0(title_prefix, " \u2013 repetition #", rep)
#               }
#               if (file_format == "png") {
#                 file_prefix <- sub("% of reads", "%%", title_prefixes[[j]], fixed = TRUE)
#                 file_prefix <- paste0(letters[[j]], ") ", file_prefix)
#                 if (!(is_all)) {
#                   file_prefix <- paste0(file_prefix, " - rep #", rep)
#                 }
#                 file_name <- paste0(file_prefix, " - ", plate_name, ".png")
#                 png(file   = file.path(sub_folder_path, file_name),
#                     width  = current_width,
#                     height = current_height,
#                     units  = "in",
#                     res    = 600
#                     )
#               }
#               use_title <- paste0(title_prefix, " (", plate_name, ")")
#               use_ccs <- paste0("ccs", ccs_numbers[[i]])
#               use_df <- subsampled_list[[j]][[rep]][[use_ccs]][[df_name]]
#               use_df <- use_df[use_df[["Plate_number"]] %in% plate_number, ]
#               PlotFunction(use_df, main_title = use_title, reorder_wells = FALSE)
#               if (file_format == "png") {
#                 dev.off()
#               }
#             }
#           }
#           if (file_format == "pdf") {
#             dev.off()
#           }
#         }
#       }
#     }
#   }
# }
#
#
#




