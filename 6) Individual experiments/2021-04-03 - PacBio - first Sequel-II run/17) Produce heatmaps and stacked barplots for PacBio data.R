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
PNGs_output_directory    <- file.path(file_output_directory, "PNGs")




# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(sql2_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(sql2_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))
load(file.path(sql2_R_objects_directory, "13) Process demultiplexed reads - with subsampling.RData"))





# Prepare constants -------------------------------------------------------

# For sparse graphics
margin_ratio <- 10
sparse_width <- 3.0 * (1 + (2 / margin_ratio))
sparse_height <- 1.5 * (2 + (3 / margin_ratio)) - (0.2 * 2)




# Export individual graphics ----------------------------------------------

ccs_numbers <- c(3, 5, 7)
accuracy_percentages <- c(99, 99.9, 99.99)

use_plate_numbers <- plates_df[["Plate_number"]][order(plates_df[["Plate_rank"]])]

for (file_format in c("png", "pdf")) {
  for (reorder_wells in c(FALSE, TRUE)) {

    order_folder <- c("original order", "re-ordered")[[as.integer(reorder_wells) + 1]]
    heatmaps_folder   <- paste0("Heatmaps - ", order_folder)
    barplots_folder   <- paste0("Stacked barplots - ", order_folder)
    sand_charts_folder <- "Sand charts"
    sparse_folder <- "Heatmaps - original - sparse"

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

        sel_name <- paste0("CCS", ccs_numbers[[i]],
                           " (", accuracy_percentages[[i]], ") - ",
                           c("i) unfiltered", "ii) filtered", "iii) filtered gRNAs")[[filter_stage]]
                           )

        if (file_format == "pdf") {
          pdf(file = file.path(plots_output_directory, heatmaps_folder, paste0("Heatmaps - ", sel_name, ".pdf")),
              width = use_width, height = use_height
              )
        } else if (file_format == "png") {
          sub_folder_path <- file.path(PNGs_output_directory, heatmaps_folder, sel_name)
          dir.create(sub_folder_path, showWarnings = FALSE)
        }
        for (plate_number in use_plate_numbers) {
          sub_df <- use_summary_df[use_summary_df[["Plate_number"]] %in% plate_number, ]
          sg_sequences_df <- library_df[library_df[["Plate_number"]] %in% plate_number, ]
          if (file_format == "png") {
            file_name <- paste0(sel_name, " - ", PlateFileName(plate_number), ".png")
            png(file   = file.path(sub_folder_path, file_name),
                width  = use_width,
                height = use_height,
                units  = "in",
                res    = 600
                )
          }
          DrawAccuracyHeatmap(sub_df,
                              main_title = titles_list[[plate_number]],
                              reorder_wells = reorder_wells
                              )
          if (file_format == "png") {
            dev.off()
          }
        }
        if (file_format == "pdf") {
          dev.off()
        }


        if (reorder_wells) {
          if (file_format == "pdf") {
            pdf(file = file.path(plots_output_directory, sand_charts_folder, paste0("Sand chart - ", sel_name, ".pdf")),
                width = 3.8, height = 6.7
                )
          } else if (file_format == "png") {
            sub_folder_path <- file.path(PNGs_output_directory, sand_charts_folder, sel_name)
            dir.create(sub_folder_path, showWarnings = FALSE)
          }
          for (plate_number in use_plate_numbers) {
            sub_df <- use_summary_df[use_summary_df[["Plate_number"]] %in% plate_number, ]
            sg_sequences_df <- library_df[library_df[["Plate_number"]] %in% plate_number, ]
            if (file_format == "png") {
              file_name <- paste0(sel_name, " - ", PlateFileName(plate_number), ".png")
              png(file   = file.path(sub_folder_path, file_name),
                  width  = 3.8,
                  height = 6.7,
                  units  = "in",
                  res    = 600
                  )
            }
            DrawReorderedSandPlots(sub_df, main_title = titles_list[[plate_number]])
            if (file_format == "png") {
              dev.off()
            }
          }
          if (file_format == "pdf") {
            dev.off()
          }
        } else {

          if (file_format == "pdf") {
            pdf(file = file.path(plots_output_directory, sparse_folder, paste0("Heatmap - ", sel_name, ".pdf")),
                width = sparse_width, height = sparse_height
                )
          } else if (file_format == "png") {
            sub_folder_path <- file.path(PNGs_output_directory, sparse_folder, sel_name)
            dir.create(sub_folder_path, showWarnings = FALSE)
          }
          for (plate_number in use_plate_numbers) {
            sub_df <- use_summary_df[use_summary_df[["Plate_number"]] %in% plate_number, ]
            sg_sequences_df <- library_df[library_df[["Plate_number"]] %in% plate_number, ]
            if (file_format == "png") {
              file_name <- paste0(sel_name, " - ", PlateFileName(plate_number), ".png")
              png(file   = file.path(sub_folder_path, file_name),
                  width  = sparse_width,
                  height = sparse_height,
                  units  = "in",
                  res    = 600
                  )
            }
            SparseAccuracyHeatmap(sub_df, main_title = titles_list[[plate_number]])
            if (file_format == "png") {
              dev.off()
            }
          }
          if (file_format == "pdf") {
            dev.off()
          }
        }

        alterations_factor <- 0.7
        alterations_width <- use_width * alterations_factor
        alterations_height <- use_height * alterations_factor
        if (file_format == "pdf") {
          pdf(file = file.path(plots_output_directory, barplots_folder, paste0("Stacked barplots - ", sel_name, ".pdf")),
              height = alterations_height, width = alterations_width
              )
        } else if (file_format == "png") {
          sub_folder_path <- file.path(PNGs_output_directory, barplots_folder, sel_name)
          dir.create(sub_folder_path, showWarnings = FALSE)
        }
        for (plate_number in use_plate_numbers) {
          sub_df <- use_summary_df[use_summary_df[["Plate_number"]] %in% plate_number, ]
          sg_sequences_df <- library_df[library_df[["Plate_number"]] %in% plate_number, ]
          if (file_format == "png") {
            file_name <- paste0(sel_name, " - ", PlateFileName(plate_number), ".png")
            png(file   = file.path(sub_folder_path, file_name),
                width  = alterations_width,
                height = alterations_height,
                units  = "in",
                res    = 600
                )
          }
          ModifiedAlterationBarplot(sub_df,
                                    main_title = titles_list[[plate_number]],
                                    reorder_wells = reorder_wells
                                    )
          if (file_format == "png") {
            dev.off()
          }
        }
        if (file_format == "pdf") {
          dev.off()
        }
      }
    }
  }
}



# Export the results of subsampling simulations ---------------------------

sampling_levels <- names(subsampled_list)
are_all <- sampling_levels == "100% sampled"
title_prefixes <- ifelse(are_all,
                         "All reads",
                         paste0("Sampled ", sub(" sampled", "", sampling_levels, fixed = TRUE), " of reads")
                         )
use_num_reps <- 3L

for (file_format in c("png", "pdf")) {
  for (i in seq_along(ccs_numbers)) {
    for (filter_stage in 1:2) {

      df_name <- c("original_summary_df", "filtered_summary_df")[[filter_stage]] # "filtered_gRNAs_df"
      use_summary_df <- use_df_list[[df_name]]
      sel_name <- paste0("CCS", ccs_numbers[[i]],
                         " (", accuracy_percentages[[i]], ") - ",
                         c("i) unfiltered", "ii) filtered", "iii) filtered gRNAs")[[filter_stage]]
                         )

      for (make_plot in c("heatmap", "sparse heatmap", "stacked barplot")) {

        if (make_plot == "heatmap") {
          folder_name <- "Heatmaps"
          PlotFunction <- DrawAccuracyHeatmap
          current_width <- use_width
          current_height <- use_height
        } else if (make_plot == "sparse heatmap") {
          folder_name <- "Heatmaps (sparse)"
          PlotFunction <- function(summary_df, main_title, reorder_wells) SparseAccuracyHeatmap(summary_df, main_title)
          current_width <- sparse_width
          current_height <- sparse_height
        } else if (make_plot == "stacked barplot") {
          folder_name <- "Stacked barplots"
          PlotFunction <- ModifiedAlterationBarplot
          current_width <- alterations_width
          current_height <- alterations_height
        }

        folder_path <- file.path("Sub-sampled", folder_name, sel_name)
        if (file_format == "pdf") {
          folder_path <- file.path(plots_output_directory, folder_path)
        } else if (file_format == "png") {
          folder_path <- file.path(PNGs_output_directory, folder_path)
        }
        dir.create(folder_path, showWarnings = FALSE)

        for (plate_number in use_plate_numbers) {
          plate_name <- paste0("Plate", formatC(plate_number, width = 2, flag = "0"),
                               " - ",
                               plates_df[["Plate_name"]][plates_df[["Plate_number"]] == plate_number]
                               )

          sg_sequences_df <- library_df[library_df[["Plate_number"]] %in% plate_number, ]

          if (file_format == "pdf") {
            pdf(file = file.path(folder_path, paste0(plate_name, ".pdf")),
                width = current_width, height = current_height
                )
          } else if (file_format == "png") {
            sub_folder_path <- file.path(folder_path, plate_name)
            dir.create(sub_folder_path, showWarnings = FALSE)
          }
          for (j in seq_along(subsampled_list)) {
            for (rep in seq_len(min(length(subsampled_list[[j]]), use_num_reps))) {
              title_prefix <- title_prefixes[[j]]
              is_all <- names(subsampled_list)[[j]] == "100% sampled"
              if (!(is_all)) {
                title_prefix <- paste0(title_prefix, " \u2013 repetition #", rep)
              }
              if (file_format == "png") {
                file_prefix <- sub("% of reads", "%%", title_prefixes[[j]], fixed = TRUE)
                file_prefix <- paste0(letters[[j]], ") ", file_prefix)
                if (!(is_all)) {
                  file_prefix <- paste0(file_prefix, " - rep #", rep)
                }
                file_name <- paste0(file_prefix, " - ", plate_name, ".png")
                png(file   = file.path(sub_folder_path, file_name),
                    width  = current_width,
                    height = current_height,
                    units  = "in",
                    res    = 600
                    )
              }
              use_title <- paste0(title_prefix, " (", plate_name, ")")
              use_ccs <- paste0("ccs", ccs_numbers[[i]])
              use_df <- subsampled_list[[j]][[rep]][[use_ccs]][[df_name]]
              use_df <- use_df[use_df[["Plate_number"]] %in% plate_number, ]
              PlotFunction(use_df, main_title = use_title, reorder_wells = FALSE)
              if (file_format == "png") {
                dev.off()
              }
            }
          }
          if (file_format == "pdf") {
            dev.off()
          }
        }
      }
    }
  }
}







