### 7th December 2021 ###



# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2r4_directory           <- file.path(experiments_directory, "2021-12-07 - fourth Sequel-II run")

s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
s2r4_R_objects_directory <- file.path(s2r4_directory, "3) R objects")



# Load data ---------------------------------------------------------------

load(file.path(s2r2_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(s2r4_R_objects_directory, "01) Process and export plate barcodes.RData"))




# Filter library_df -------------------------------------------------------

are_present <- library_df[["Plate_name"]] %in% plates_df[["Plate_name"]]
library_df <- library_df[are_present, ]
row.names(library_df) <- NULL




# Save data ---------------------------------------------------------------

save(list = "library_df",
     file = file.path(s2r4_R_objects_directory, "03) Import and process sgRNA sequences.RData")
     )



