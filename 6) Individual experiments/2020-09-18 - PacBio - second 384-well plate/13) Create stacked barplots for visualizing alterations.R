### 21st October 2020 ###




# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "09) Producing heatmaps.R"))
source(file.path(R_functions_directory, "11) Creating stacked barplots for visualizing alterations.R"))




# Define folder paths -----------------------------------------------------

plate2_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-18 - PacBio - second 384-well plate")
p2_R_objects_directory <- file.path(plate2_directory, "2) R objects")
file_output_directory  <- file.path(plate2_directory, "3) Output")
plots_output_directory <- file.path(file_output_directory, "Figures")




# Load data ---------------------------------------------------------------

load(file.path(p2_R_objects_directory, "01) Import and process sgRNA sequences.RData"))
load(file.path(p2_R_objects_directory, "07) Process demultiplexed PacBio reads.RData"))





# Exclude 4 problematic wells, for the time being -------------------------

sg_sequences_df[["Empty_well"]] <- ifelse(sg_sequences_df[["Well_number"]] == 2,
                                          TRUE, FALSE
                                          )




# Set up loop -------------------------------------------------------------

for (smrtlink_version in c(7, 9)) {

  for (reorder_wells in c(FALSE, TRUE)) {

    version_folder <- paste0("SmrtLink ", smrtlink_version)
    order_folder <- c("original order", "re-ordered")[[as.integer(reorder_wells) + 1]]
    order_folder <- paste0("Stacked barplots - ", order_folder)
    plots_dir <- file.path(plots_output_directory, version_folder, order_folder)


    # Draw the accuracy plots in the console ----------------------------------

    DrawAlterationBarplot(sl7_ccs5_df_list[["original_summary_df"]],
                          main_title = ccs5_title, reorder_wells = reorder_wells
                          )
    DrawAlterationBarplot(sl7_ccs5_df_list[["filtered_summary_df"]],
                          main_title = ccs5_title, reorder_wells = reorder_wells
                          )
    DrawAlterationBarplot(sl7_ccs5_df_list[["filtered_gRNAs_df"]],
                          main_title = ccs5_title, reorder_wells = reorder_wells
                          )

    DrawAlterationBarplot(sl7_ccs3_df_list[["original_summary_df"]],
                          main_title = ccs3_title, reorder_wells = reorder_wells
                          )
    DrawAlterationBarplot(sl7_ccs3_df_list[["filtered_summary_df"]],
                          main_title = ccs3_title, reorder_wells = reorder_wells
                          )
    DrawAlterationBarplot(sl7_ccs3_df_list[["filtered_gRNAs_df"]],
                          main_title = ccs3_title, reorder_wells = reorder_wells
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
      DrawAlterationBarplot(summary_df, main_title = main_title,
                            reorder_wells = reorder_wells
                            )
      dev.off()
    }

    file_name_prefix <- paste0("Alteration barplots - SmrtLink ", smrtlink_version)

    SavePNG(sl7_ccs5_df_list[["original_summary_df"]],
            paste0(file_name_prefix, " - CCS5 (99.9) - original"),
            main_title = ccs5_title
            )
    SavePNG(sl7_ccs5_df_list[["filtered_summary_df"]],
            paste0(file_name_prefix, " - CCS5 (99.9) - filtered"),
            main_title = ccs5_title
            )
    SavePNG(sl7_ccs5_df_list[["filtered_gRNAs_df"]],
            paste0(file_name_prefix, " - CCS5 (99.9) - filtered gRNAs"),
            main_title = ccs5_title
            )

    SavePNG(sl7_ccs3_df_list[["original_summary_df"]],
            paste0(file_name_prefix, " - CCS3 (99) - original"),
            main_title = ccs3_title
            )
    SavePNG(sl7_ccs3_df_list[["filtered_summary_df"]],
            paste0(file_name_prefix, " - CCS3 (99) - filtered"),
            main_title = ccs3_title
            )
    SavePNG(sl7_ccs3_df_list[["filtered_gRNAs_df"]],
            paste0(file_name_prefix, " - CCS3 (99) - filtered gRNAs"),
            main_title = ccs3_title
            )



    # Produce the accuracy PDF ------------------------------------------------

    pdf(file = file.path(plots_dir, paste0(file_name_prefix, " - original.pdf")),
        height = use_height,
        width  = use_width
        )
    DrawAlterationBarplot(sl7_ccs5_df_list[["original_summary_df"]],
                          main_title = ccs5_title, reorder_wells = reorder_wells
                          )
    DrawAlterationBarplot(sl7_ccs3_df_list[["original_summary_df"]],
                          main_title = ccs3_title, reorder_wells = reorder_wells
                          )
    dev.off()


    pdf(file = file.path(plots_dir, paste0(file_name_prefix, " - filtered.pdf")),
        height = use_height,
        width  = use_width
        )
    DrawAlterationBarplot(sl7_ccs5_df_list[["filtered_summary_df"]],
                          main_title = ccs5_title, reorder_wells = reorder_wells
                          )
    DrawAlterationBarplot(sl7_ccs3_df_list[["filtered_summary_df"]],
                          main_title = ccs3_title, reorder_wells = reorder_wells
                          )
    dev.off()


    pdf(file = file.path(plots_dir, paste0(file_name_prefix, " - filtered gRNAs.pdf")),
        height = use_height,
        width  = use_width
        )
    DrawAlterationBarplot(sl7_ccs5_df_list[["filtered_gRNAs_df"]],
                          main_title = ccs5_title, reorder_wells = reorder_wells
                          )
    DrawAlterationBarplot(sl7_ccs3_df_list[["filtered_gRNAs_df"]],
                          main_title = ccs3_title, reorder_wells = reorder_wells
                          )
    dev.off()


    # End loop ----------------------------------------------------------------
  }
}


