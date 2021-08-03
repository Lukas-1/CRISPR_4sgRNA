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
PNGs_output_directory    <- file.path(file_output_directory, "PNGs", "Schematics")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(sql2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(sql2_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(sql2_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))





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

use_plate_numbers <- plates_df[["Plate_number"]][order(plates_df[["Plate_rank"]])]

DrawSchematicsForAllPlates()














