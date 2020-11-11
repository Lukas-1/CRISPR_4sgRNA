### 28 May 2020 ###




# Import packages and source code -----------------------------------------

CRISPR_root_directory  <- "~/CRISPR"
file_directory         <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory  <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "09) Producing heatmaps.R"))




# Define folder paths -----------------------------------------------------

R_objects_directory    <- file.path(file_directory, "3) R objects")
file_output_directory  <- file.path(file_directory, "5) Output")
plots_output_directory <- file.path(file_output_directory, "Figures")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(R_objects_directory, "09) Process demultiplexed PacBio reads.RData"))





# Set up loop -------------------------------------------------------------

for (smrtlink_version in c(7, 9)) {

  for (reorder_wells in c(FALSE, TRUE)) {

    version_folder <- paste0("SmrtLink ", smrtlink_version)
    order_folder <- c("original order", "re-ordered")[[as.integer(reorder_wells) + 1]]
    order_folder <- paste0("Heatmaps - ", order_folder)
    plots_dir <- file.path(plots_output_directory, version_folder, order_folder)

    if (smrtlink_version == 7) {
      ccs3_df_list <- sl7_ccs3_df_list
      ccs5_df_list <- sl7_ccs5_df_list
    } else {
      ccs3_df_list <- sl9_ccs3_df_list
      ccs5_df_list <- sl9_ccs5_df_list
    }

    # Draw the accuracy plots in the console ----------------------------------

    DrawAccuracyHeatmap(ccs5_df_list[["original_summary_df"]],
                        main_title = ccs5_title, reorder_wells = reorder_wells
                        )
    DrawAccuracyHeatmap(ccs5_df_list[["filtered_summary_df"]],
                        main_title = ccs5_title, reorder_wells = reorder_wells
                        )
    DrawAccuracyHeatmap(ccs5_df_list[["filtered_gRNAs_df"]],
                        main_title = ccs5_title, reorder_wells = reorder_wells
                        )

    DrawAccuracyHeatmap(sl7_ccs3_df_list[["original_summary_df"]],
                        main_title = ccs3_title, reorder_wells = reorder_wells
                        )
    DrawAccuracyHeatmap(ccs3_df_list[["filtered_summary_df"]],
                        main_title = ccs3_title, reorder_wells = reorder_wells
                        )
    DrawAccuracyHeatmap(ccs3_df_list[["filtered_gRNAs_df"]],
                        main_title = ccs5_title, reorder_wells = reorder_wells
                        )


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

    SavePNG(ccs5_df_list[["original_summary_df"]],
            paste0(file_name_prefix, " - CCS5 (99.9) - original"),
            main_title = ccs5_title
            )
    SavePNG(ccs5_df_list[["filtered_summary_df"]],
            paste0(file_name_prefix, " - CCS5 (99.9) - filtered"),
            main_title = ccs5_title
            )
    SavePNG(ccs5_df_list[["filtered_gRNAs_df"]],
            paste0(file_name_prefix, " - CCS5 (99.9) - filtered gRNAs"),
            main_title = ccs5_title
            )

    SavePNG(ccs3_df_list[["original_summary_df"]],
            paste0(file_name_prefix, " - CCS3 (99) - original"),
            main_title = ccs3_title
            )
    SavePNG(ccs3_df_list[["filtered_summary_df"]],
            paste0(file_name_prefix, " - CCS3 (99) - filtered"),
            main_title = ccs3_title
            )
    SavePNG(ccs3_df_list[["filtered_gRNAs_df"]],
            paste0(file_name_prefix, " - CCS3 (99) - filtered gRNAs"),
            main_title = ccs3_title
            )



    # Produce the accuracy PDF ------------------------------------------------

    pdf(file = file.path(plots_dir, paste0(file_name_prefix, " - original.pdf")),
        height = use_height,
        width  = use_width
        )
    DrawAccuracyHeatmap(ccs5_df_list[["original_summary_df"]],
                        main_title = ccs5_title, reorder_wells = reorder_wells
                        )
    DrawAccuracyHeatmap(ccs3_df_list[["original_summary_df"]],
                        main_title = ccs3_title, reorder_wells = reorder_wells
                        )
    dev.off()


    pdf(file = file.path(plots_dir, paste0(file_name_prefix, " - filtered.pdf")),
        height = use_height,
        width  = use_width
        )
    DrawAccuracyHeatmap(ccs5_df_list[["filtered_summary_df"]],
                        main_title = ccs5_title, reorder_wells = reorder_wells
                        )
    DrawAccuracyHeatmap(ccs3_df_list[["filtered_summary_df"]],
                        main_title = ccs3_title, reorder_wells = reorder_wells
                        )

    dev.off()


    pdf(file = file.path(plots_dir, paste0(file_name_prefix, " - filtered gRNAs.pdf")),
        height = use_height,
        width  = use_width
        )
    DrawAccuracyHeatmap(ccs5_df_list[["filtered_gRNAs_df"]],
                        main_title = ccs5_title, reorder_wells = reorder_wells
                        )
    DrawAccuracyHeatmap(ccs3_df_list[["filtered_gRNAs_df"]],
                        main_title = ccs3_title, reorder_wells = reorder_wells
                        )

    dev.off()


    # End loop ----------------------------------------------------------------
  }
}



