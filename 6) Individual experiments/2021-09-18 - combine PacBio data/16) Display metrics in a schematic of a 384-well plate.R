### 30th September 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "01) Define titles and labels.R"))
source(file.path(R_functions_directory, "16) Showing metrics in a schematic of a 384-well plate.R"))




# Define folder paths -----------------------------------------------------

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2rC_directory           <- file.path(experiments_directory, "2021-09-18 - combine PacBio data")
p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")
s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
s2rC_R_objects_directory <- file.path(s2rC_directory, "3) R objects")
file_output_directory    <- file.path(s2rC_directory, "5) Output")
plots_output_directory   <- file.path(file_output_directory, "Figures", "Schematics of a 384-well plate")
# PNGs_output_directory    <- file.path(file_output_directory, "PNGs", "Schematics")




# Load data ---------------------------------------------------------------

load(file.path(s2r2_R_objects_directory, "30) Calculate correction factors to account for mean read counts.RData"))

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(s2rC_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2rC_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(s2rC_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))




# General preparation -----------------------------------------------------

sg_sequences_df[["Empty_well"]] <- FALSE

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




# Export all plates -------------------------------------------------------

use_plate_numbers <- plates_df[["Plate_number"]]#[order(plates_df[["Plate_rank"]])]

DrawSchematicsForAllPlates(export_PNGs = FALSE)






