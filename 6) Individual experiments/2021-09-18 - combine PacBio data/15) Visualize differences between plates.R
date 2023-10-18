### 30th September 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory      <- "~/CRISPR_4sgRNA"
experiments_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory           <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
s2r1_directory             <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
p1_R_functions_directory   <- file.path(plate1_directory, "1) R functions")
s2r1_R_functions_directory <- file.path(s2r1_directory, "1) R functions")

source(file.path(p1_R_functions_directory, "01) Define titles and labels.R"))
source(file.path(s2r1_R_functions_directory, "04) Visualizing differences between plates.R"))




# Define folder paths -----------------------------------------------------

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2rC_directory           <- file.path(experiments_directory, "2021-09-18 - combine PacBio data")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
s2rC_R_objects_directory <- file.path(s2rC_directory, "3) R objects")
file_output_directory    <- file.path(s2rC_directory, "5) Output")
plots_output_directory   <- file.path(file_output_directory, "Figures", "Compare plates")




# Load data ---------------------------------------------------------------

load(file.path(s2r2_R_objects_directory, "30) Calculate correction factors to account for mean read counts.RData"))

load(file.path(s2rC_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2rC_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))




# Make preparations -------------------------------------------------------

titles_list <- c(list("Count_total" = "Number of reads per well"),
                 titles_list
                 )

plate_numbers <- unique(ccs5_df_list[["original_summary_df"]][["Plate_number"]])
plates_df <- plates_df[plates_df[["Plate_number"]] %in% plate_numbers, ]
row.names(plates_df) <- NULL

matches_vec <- match(plates_df[, "Plate_name"], extended_df[, "Plate_name"])
pools_vec <- extended_df[matches_vec, "Run3_pool"]
plates_df[, "Highlight_color"] <- ifelse(pools_vec %in% 3,
                                         brewer.pal(9, "Blues")[[7]],
                                         ifelse(pools_vec %in% 4,
                                                brewer.pal(9, "Purples")[[7]],
                                                "#000000"
                                                )
                                         )




# Draw example plots ------------------------------------------------------

ComparePlates(ccs7_df_list[["filtered_summary_df"]], "Count_total",
              use_cex = 0.075, side_space = -2
              )




# Export graphics ---------------------------------------------------------

DrawAllPlateComparisons(use_cex          = 0.175,
                        beeswarm_spacing = 0.3,
                        beeswarm_corral  = "omit",
                        side_space       = -3.5,
                        use_width        = 22,
                        use_height       = 6.5,
                        export_PNGs      = FALSE,
                        exclude_controls = TRUE
                        )




