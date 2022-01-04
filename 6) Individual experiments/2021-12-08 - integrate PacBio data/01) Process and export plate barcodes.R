### 8th December 2021 ###



# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")

s2rC_directory           <- file.path(experiments_directory, "2021-09-18 - combine PacBio data")
s2rI_directory           <- file.path(experiments_directory, "2021-12-08 - integrate PacBio data")

s2rC_R_objects_directory <- file.path(s2rC_directory, "3) R objects")
s2rI_R_objects_directory <- file.path(s2rI_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(s2rC_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2rC_R_objects_directory, "03) Import and process sgRNA sequences.RData"))




# Create the new plates_df ------------------------------------------------

exclude_columns <- c("Run2_pool", "Number_96wp")
plates_df <- plates_df[, !(names(plates_df) %in% exclude_columns)]
plates_df <- plates_df[plates_df[, "Plate_number"] %in% library_df[, "Plate_number"], ]
row.names(plates_df) <- NULL




# Save data ---------------------------------------------------------------

save(list = "plates_df",
     file = file.path(s2rI_R_objects_directory, "01) Process and export plate barcodes.RData")
     )



