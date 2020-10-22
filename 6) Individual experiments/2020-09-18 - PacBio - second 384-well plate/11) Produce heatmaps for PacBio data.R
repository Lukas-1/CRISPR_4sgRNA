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




# Load data ---------------------------------------------------------------

load(file.path(p2_R_objects_directory, "01) Import and process sgRNA sequences.RData"))
load(file.path(p2_R_objects_directory, "08) Process demultiplexed PacBio reads.RData"))




# Exclude 4 problematic wells, for the time being -------------------------

sg_sequences_df[["Empty_well"]] <- ifelse(sg_sequences_df[["Well_number"]] %in% c(2, 171, 284, 285),
                                         TRUE, FALSE
                                         )




# Define plot titles ------------------------------------------------------

ccs3_title <- expression(plain({"Long-read sequencing of plasmids" *
                               " (" >= "3 consensus reads "} *
                               "and " >= "99% accuracy)"
                               ))

ccs5_title <- expression(plain({"Long-read sequencing of plasmids" *
                               " (" >= "5 consensus reads "} *
                               "and " >= "99.9% accuracy)"
                               ))



# Define plot dimensions --------------------------------------------------

use_height <- 7
use_width <- 6.5




# Set up loop -------------------------------------------------------------

for (smrtlink_version in c(7)) {
# for (smrtlink_version in c(7, 9)) {

  for (reorder_wells in c(FALSE, TRUE)) {

    version_folder <- paste0("SmrtLink ", smrtlink_version)
    order_folder <- c("original order", "re-ordered")[[as.integer(reorder_wells) + 1]]
    order_folder <- paste0("Heatmaps - ", order_folder)
    plots_dir <- file.path(plots_output_directory, version_folder, order_folder)


    # Draw the accuracy plots in the console ----------------------------------

    DrawAccuracyHeatmap(sl7_ccs5_df_list[["original_summary_df"]],
                        main_title = ccs5_title, reorder_wells = reorder_wells
                        )
    DrawAccuracyHeatmap(sl7_ccs5_df_list[["filtered_summary_df"]],
                        main_title = ccs5_title, reorder_wells = reorder_wells
                        )

    # DrawAccuracyHeatmap(sl7_ccs3_df_list[["original_summary_df"]],
    #                     main_title = ccs3_title, reorder_wells = reorder_wells
    #                     )
    # DrawAccuracyHeatmap(sl7_ccs3_df_list[["filtered_summary_df"]],
    #                     main_title = ccs3_title, reorder_wells = reorder_wells
    #                     )



    # Produce the accuracy PNGs -----------------------------------------------

    SavePNG <- function(summary_df, file_name, main_title) {
      full_path <- file.path(plots_dir, paste0(file_name, ".png"))
      png(filename = full_path,
          res      = 600,
          height   = use_height,
          width    = use_width,
          units    = "in"
      )
      DrawAccuracyHeatmap(summary_df, main_title = main_title,
                          reorder_wells = reorder_wells
                          )
      dev.off()
    }

    file_name_prefix <- paste0("Accuracy heatmap - SmrtLink ", smrtlink_version)

    SavePNG(sl7_ccs5_df_list[["original_summary_df"]],
            paste0(file_name_prefix, " - CCS5 (99.9) - original"),
            main_title = ccs5_title
            )
    SavePNG(sl7_ccs5_df_list[["filtered_summary_df"]],
            paste0(file_name_prefix, " - CCS5 (99.9) - filtered"),
            main_title = ccs5_title
            )

    # SavePNG(sl7_ccs3_df_list[["original_summary_df"]],
    #         paste0(file_name_prefix, " - CCS3 (99) - original"),
    #         main_title = ccs3_title
    #         )
    # SavePNG(sl7_ccs3_df_list[["filtered_summary_df"]],
    #         paste0(file_name_prefix, " - CCS3 (99) - filtered"),
    #         main_title = ccs3_title
    #         )




    # Produce the accuracy PDF ------------------------------------------------

    pdf(file = file.path(plots_dir, paste0(file_name_prefix, " - original.pdf")),
        height = use_height,
        width  = use_width
        )
    DrawAccuracyHeatmap(sl7_ccs5_df_list[["original_summary_df"]],
                        main_title = ccs5_title, reorder_wells = reorder_wells
                        )
    # DrawAccuracyHeatmap(sl7_ccs3_df_list[["original_summary_df"]],
    #                     main_title = ccs3_title, reorder_wells = reorder_wells
    #                     )
    dev.off()


    pdf(file = file.path(plots_dir, paste0(file_name_prefix, " - filtered.pdf")),
        height = use_height,
        width  = use_width
        )
    DrawAccuracyHeatmap(sl7_ccs5_df_list[["filtered_summary_df"]],
                        main_title = ccs5_title, reorder_wells = reorder_wells
                        )
    # DrawAccuracyHeatmap(sl7_ccs3_df_list[["filtered_summary_df"]],
    #                     main_title = ccs3_title, reorder_wells = reorder_wells
    #                     )
    dev.off()



    # End loop ----------------------------------------------------------------
  }
}



