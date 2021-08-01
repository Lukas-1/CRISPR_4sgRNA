### 29th July 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "01) Define titles and labels.R"))
source(file.path(R_functions_directory, "16) Showing metrics in a schematic of a 384-well plate.R"))




# Define folder paths -----------------------------------------------------

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
file_output_directory    <- file.path(s2r2_directory, "5) Output")
plots_output_directory   <- file.path(file_output_directory, "Figures", "Schematics of a 384-well plate")
PNGs_output_directory    <- file.path(file_output_directory, "PNGs", "Schematics")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(s2r2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2r2_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(s2r2_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))





# General preparation -----------------------------------------------------

sg_sequences_df[["Empty_well"]] <- FALSE





# Display metrics in the layout of a 384-well plate -----------------------

use_df <- ccs7_df_list[["filtered_summary_df"]]

this_plate <- 24L

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

BarPlotPanel(summary_sub_df,
             "Count_all_4",
             sg_sub_df,
             show_low_read_numbers = TRUE,
             outline_few_reads = TRUE
             )

BarPlotPanel(summary_sub_df,
             "Count_pr_all_4",
             sg_sub_df,
             show_low_read_numbers = TRUE,
             outline_few_reads = TRUE
             )


# Export all plates -------------------------------------------------------

use_plate_numbers <- plates_df[["Plate_number"]]#[order(plates_df[["Plate_rank"]])]

ccs_numbers <- c(3, 5, 7)
accuracy_percentages <- c(99, 99.9, 99.99)

count_metrics <- c(
  # "Count_sg1_cr1", "Count_sg2_cr2", "Count_sg3_cr3", "Count_sg4_cr4",
  "Count_mean_sg1to4",
  "Count_at_least_1", "Count_at_least_2", "Count_at_least_3", "Count_all_4",
  "Count_pr_at_least_1", "Count_pr_at_least_2", "Count_pr_at_least_3", "Count_pr_all_4",
  "Count_all_4_promoters", "Count_whole_plasmid"
)

percentages_metrics <- c(
  "Num_contaminated_reads", "Num_under_2kb",
  count_metrics,
  "Num_reads_with_deletions_exceeding_20bp",
  "Num_reads_with_sgRNA_deletion",
  "Num_reads_with_deletions_spanning_tracrRNAs",
  "Num_reads_with_deletions_spanning_promoters",
  "Num_cross_plate_contaminated",
  "Num_contaminated_reads_aligned"
)

binarized_metrics <- c(
  "Binary_all_four_guides",
  paste0("Binary_", tolower(substr(count_metrics, 1, 1)),
         substr(count_metrics, 2, nchar(count_metrics))
         ),
  "Binary_num_contaminated_reads",
  "Binary_num_reads_with_deletions_exceeding_20bp",
  "Binary_num_reads_with_deletions_spanning_tracrRNAs"
)


plate_labels <- paste0("Plate #", plates_df[["Plate_number"]], " \u2013 ", plates_df[["Plate_name"]])
plate_number_width <- max(nchar(as.character(plates_df[["Plate_number"]])))

# for (file_format in c("png", "pdf")) {
for (file_format in c("pdf")) {
  for (i in seq_along(ccs_numbers)) {
    use_df_list <- get(paste0("ccs", ccs_numbers[[i]], "_df_list"))
    for (filter_stage in 1:2) {
      df_name <- c("original_summary_df", "filtered_summary_df")[[filter_stage]] # "filtered_gRNAs_df"

      use_summary_df <- use_df_list[[df_name]]

      folder_name <- paste0("CCS", ccs_numbers[[i]],
                            " (", accuracy_percentages[[i]], ") - ",
                            c("i) unfiltered", "ii) filtered", "iii) filtered gRNAs")[[filter_stage]]
                            )
      if (file_format == "png") {
        folder_path <- file.path(PNGs_output_directory, folder_name)
      } else {
        folder_path <- file.path(plots_output_directory, folder_name)
      }

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

            current_metric_name <- sub("^Num_reads_with_", "", current_metrics[[j]])
            current_metric_name <- sub("^(Num|Count)_", "", current_metric_name)
            current_metric_name <- sub("^Binary_num_reads_with_", "Binary_", current_metric_name)
            current_metric_name <- sub("^Binary_(num|count)_", "Binary_", current_metric_name)
            current_metric_name <- paste0(formatC(j, width = 2, flag = "0"),
                                          ") ",  current_metric_name
                                          )

            expansion_factor <- 5
            use_width <- 2 * expansion_factor
            use_height <- 1 * expansion_factor

            if (file_format == "pdf") {
              pdf_name <- paste0(file_prefix, " - ", current_metric_name, ".pdf")
              pdf(file   = file.path(sub_folder_path, pdf_name),
                  width  = use_width,
                  height = use_height
                  )
            } else if (file_format == "png") {
              metric_path <- file.path(sub_folder_path, current_metric_name)
              dir.create(metric_path, showWarnings = FALSE)
            }

            for (plate_number in use_plate_numbers) {
              if (file_format == "png") {
                plate_name <- paste0(formatC(plate_number, width = plate_number_width, flag = "0"), ") - ",
                                     plates_df[["Plate_name"]][plates_df[["Plate_number"]] == plate_number]
                                     )
                file_name <- paste0(file_prefix, " - ", plate_name, ".png")
                png(file   = file.path(metric_path, file_name),
                    width  = use_width,
                    height = use_height,
                    units  = "in",
                    res    = 600
                    )
              }
              summary_sub_df <- use_summary_df[use_summary_df[["Plate_number"]] %in% plate_number, ]
              sg_sub_df <- sg_sequences_df[sg_sequences_df[["Plate_number"]] %in% plate_number, ]
              BarPlotPanel(summary_sub_df,
                           current_metrics[[j]],
                           sg_sub_df,
                           number_wells          = number_wells,
                           top_text              = plate_labels[plates_df[["Plate_number"]] == plate_number],
                           show_low_read_numbers = TRUE,
                           outline_few_reads     = binarize
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
  }
}







