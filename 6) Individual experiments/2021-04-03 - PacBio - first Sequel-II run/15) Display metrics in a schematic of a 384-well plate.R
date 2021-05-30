### 19th September 2020 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "01) Define titles and labels.R"))
source(file.path(R_functions_directory, "16) Showing metrics in a schematic of a 384-well plate.R"))




# Define folder paths -----------------------------------------------------

sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")
file_output_directory    <- file.path(sql2_directory, "5) Output")
plots_output_directory   <- file.path(file_output_directory, "Figures", "Schematics of a 384-well plate")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(sql2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(sql2_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(sql2_R_objects_directory, "09) Process demultiplexed PacBio reads.RData"))





# Display metrics in the layout of a 384-well plate -----------------------

use_df <- ccs7_df_list[["filtered_summary_df"]]

this_plate <- 3

summary_sub_df <- use_df[use_df[["Plate_number"]] %in% this_plate, ]
sg_sub_df <- sg_sequences_df[sg_sequences_df[["Plate_number"]] %in% this_plate, ]

BarPlotPanel(summary_sub_df,
             "Count_at_least_1",
             sg_sub_df,
             show_low_read_numbers = TRUE
             )

BarPlotPanel(summary_sub_df,
             "Count_at_least_1",
             sg_sub_df,
             show_low_read_numbers = TRUE,
             outline_few_reads = TRUE
             )


BarPlotPanel(summary_sub_df,
             "Binary_count_at_least_1",
             sg_sub_df,
             show_low_read_numbers = TRUE,
             outline_few_reads = TRUE
             )

BarPlotPanel(summary_sub_df,
             "Binary_count_all_4",
             sg_sub_df,
             show_low_read_numbers = TRUE,
             outline_few_reads = TRUE
             )

BarPlotPanel(summary_sub_df,
             "Binary_all_four_guides",
             sg_sub_df,
             show_low_read_numbers = TRUE,
             outline_few_reads = TRUE
             )

BarPlotPanel(summary_sub_df,
             "Binary_count_mean_sg1to4",
             sg_sub_df,
             show_low_read_numbers = TRUE,
             outline_few_reads = TRUE
             )




# Export all plates -------------------------------------------------------

use_plate_numbers <- plates_df[["Plate_number"]][order(plates_df[["Plate_rank"]])]

ccs_numbers <- c(3, 5, 7)
accuracy_percentages <- c(99, 99.9, 99.99)

count_metrics <- c(
  # "Count_sg1_cr1", "Count_sg2_cr2", "Count_sg3_cr3", "Count_sg4_cr4",
  "Count_mean_sg1to4",
  "Count_at_least_1", "Count_at_least_2", "Count_at_least_3", "Count_all_4",
  "Count_all_4_promoters", "Count_whole_plasmid"
)

percentages_metrics <- c(
  "Num_contaminated_reads", "Num_under_2kb", count_metrics
)

binarized_metrics <- c(
  "Binary_all_four_guides",
  paste0("Binary_", tolower(substr(count_metrics, 1, 1)),
         substr(count_metrics, 2, nchar(count_metrics))
         )
)


for (i in seq_along(ccs_numbers)) {
  use_df_list <- get(paste0("ccs", ccs_numbers[[i]], "_df_list"))
  for (filter_stage in 1:2) {
    df_name <- c("original_summary_df", "filtered_summary_df")[[filter_stage]] # "filtered_gRNAs_df"

    use_summary_df <- use_df_list[[df_name]]

    folder_name <- paste0("CCS", ccs_numbers[[i]],
                          " (", accuracy_percentages[[i]], ") - ",
                          c("i) unfiltered", "ii) filtered", "iii) filtered gRNAs")[[filter_stage]]
                          )
    folder_path <- file.path(plots_output_directory, folder_name)
    dir.create(folder_path, showWarnings = FALSE)

    file_prefix <- folder_name

    for (number_wells in c(TRUE, FALSE)) {
      for (binarize in c(TRUE, FALSE)) {

        if (binarize) {
          current_metrics <- binarized_metrics
        } else {
          current_metrics <- percentages_metrics
        }

        for (j in seq_along(current_metrics)) {

          if (binarize) {
            sub_folder_path <- file.path(folder_path, "Binarized - ")
          } else {
            sub_folder_path <- file.path(folder_path, "Percentages - ")
          }
          if (number_wells) {
            sub_folder_path <- paste0(sub_folder_path, "numbered")
          } else {
            sub_folder_path <- paste0(sub_folder_path, "plain")
          }
          dir.create(sub_folder_path, showWarnings = FALSE)

          expansion_factor <- 5
          file_name <- paste0(file_prefix, " - ",
                              formatC(j, width = 2, flag = "0"), ") ",
                              current_metrics[[j]],
                              ".pdf"
                              )
          pdf(file   = file.path(sub_folder_path, file_name),
              width  = 2 * expansion_factor,
              height = 1 * expansion_factor
              )

          plate_labels <- paste0("Plate #", plates_df[["Plate_number"]], " \u2013 ", plates_df[["Plate_name"]])
          for (plate_number in use_plate_numbers) {
            summary_sub_df <- use_summary_df[use_summary_df[["Plate_number"]] %in% plate_number, ]
            sg_sub_df <- sg_sequences_df[sg_sequences_df[["Plate_number"]] %in% plate_number, ]
            BarPlotPanel(summary_sub_df,
                         current_metrics[[j]],
                         sg_sub_df,
                         number_wells = number_wells,
                         top_text = plate_labels[plates_df[["Plate_number"]] == plate_number],
                         show_low_read_numbers = TRUE,
                         outline_few_reads = binarize
                         )
          }
          dev.off()
        }

      }

    }
  }
}








