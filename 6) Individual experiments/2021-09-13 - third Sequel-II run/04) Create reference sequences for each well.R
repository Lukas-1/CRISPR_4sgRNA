### 14th September 2021 ###



# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")

s2r2_directory           <- file.path(experiments_directory, "2021-07-24 - second Sequel-II run")
s2r3_directory           <- file.path(experiments_directory, "2021-09-13 - third Sequel-II run")

s2r2_R_objects_directory <- file.path(s2r2_directory, "3) R objects")
s2r3_R_objects_directory <- file.path(s2r3_directory, "3) R objects")



# Load data ---------------------------------------------------------------

load(file.path(s2r2_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(s2r3_R_objects_directory, "01) Process and export plate barcodes.RData"))




# Filter sg_sequences_df --------------------------------------------------

are_present <- sg_sequences_df[["Plate_name"]] %in% plates_df[["Plate_name"]]
sg_sequences_df <- sg_sequences_df[are_present, ]
row.names(sg_sequences_df) <- NULL




# Save data ---------------------------------------------------------------

save(list = "sg_sequences_df",
     file = file.path(s2r3_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData")
     )


